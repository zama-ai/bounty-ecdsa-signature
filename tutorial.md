# Tutorial

# Intro

This tutorial provides an overview of ECDSA fundamentals, delves into the development of its homomorphic variant, discusses finite field and elliptic curve operations, and offers guidance on optimizing performance based on THFE characteristics.

# ECDSA

ECDSA, or Elliptic Curve Digital Signature Algorithm, is a widely used cryptographic algorithm. It provides a method for verifying the authenticity and integrity of digital data, such as messages or transactions. ECDSA relies on the mathematics of elliptic curves to generate digital signatures, which can be used to confirm that a message has been signed by a specific private key and has not been tampered with during transmission. It's commonly employed in secure communication and digital authentication systems, such as in cryptocurrencies like Ethereum.

This algorithm work on an elliptic curve and there’re 3 main operation for the algorithm namely

- Key Generation - KeyGen(randomness) → (public key, private key), normally this
- Signing - Sign(private key, msg, nonce) → signature
- Recovering - Recover (msg, signature) → (signer) public key
- Verifying - Verify (msg, signature, public key) → result (True, False)

In this work we will focus on the signing part. 

### Eliptic Curve Arithmetic

Generically elliptic curve will be written in the **Weierstrass** form/equation

$$
y^2=x^3+ax+b
$$

But some time we will see the curve in the other form too. For example, **Edwards**, **Jacobian**, **Montgomery**. Which we will not cover in this topic since SECP256k1 uses Weierstrass equation.

- $a, b$ - $a$ in $ax$ and $b$ in $+b$ in the curve equation which define it shape (and properties)
- $p$ - Prime Modulus, a prime number that defines the finite field over which the elliptic curve.
- $G$ - Generator point, a point $(x,y)$ on the curve that used to generate all others points.
- $n$ - Order of subgroup, aka how many number of point do we have in generator until such generator generate the generator point again. In a sense, $n$ defines range of possible private key. This number can be algorithmically found by [Schoof algorithm](https://en.wikipedia.org/wiki/Schoof%27s_algorithm#The_algorithm).

### Secp256k1

Secp256k1 is an elliptic curve that wieldy adopted for it performance (ex. Bitcoin). the parameter and good property from the following parameter 

- a = 0, b = 7 - **a = 0** eliminate the need to compute anything related to **ax** term.
- p = `FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFE FFFFFC2F` or $2^{256} - 2^{32} - 2^9 - 2^8 - 2^7 - 2^6 - 2^4 - 1$. This prime basically 2^256 minus some relatively small number. apart from being close to 256 bits which is the amount of bit that we want. It also allow some speed up trick that we’ll go into more detail in later section of this tutorial
- n = `FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFE BAAEDCE6 AF48A03B BFD25E8C D0364141`, large subgroup order.

### Signing

We’ll assume that you already have a private key and a message. 

- First hash your message $z = hash(m)$ and mod (trim) to the size of n. but in this case the size already correct so there’s no need to trim.
- Pick a secret random $k$ where $k$ is between $1$ and $n-1$
- Calculate $kG=(x,y)$ where $G$ is the prime subgroup generator of Secp256k1. Note that this operation mean scalar multiply the generator point with $k$ and all of the inside operation is calculated using $\mod p$.
- Find $r = x \mod n$. if $r = 0$ then back to pick $k$.
- Find $s = k^{-1}\cdot(z + r \cdot sk) \mod n$ where $sk$ is the private key of the user, and if $s = 0$ then back to pick $k$. Note that $k^{-1}$ is not simple arithmetic inverse $1\div{k}$ but is a **modular inverse** where $k\cdot k^{-1}=1\mod n$.
- $(r,s)$ is the signature

# Implementation

To implement homomorphic ECDSA signing, we will begin by working with small building blocks in the most basic implementation. From there, we will gradually and systematically expand upon these foundations, adding complexity to our project.

There’s three main category of operation involving the ECDDSA algorithm.

- Primitive homomorphic operation (exposed by `tfhe-rs`)
- Finite field operation
- Elliptic curve operation, which depends on finite field operation.

## Primitive Homomorphic Operation

`tfhe-rs` exposes many primitive operators for us to use. Such as `add`, `sub`, `mul`, `div`, `gt`, `lt`, `eq`, etc. Which the computation cost and time are different from their normal variant. So, we need to briefly know which operator to prioritize if we want the algorithm to run fastest.

We benchmark each operation and came in conclusion that (in terms of speed)

`All bits operation` << `comparators` << `add, sub` <<<< `mul` <<<< `div`

Note that computation time also grows in exponential to bit length, and trivial or plain text operand are faster than cypher text.

```rust
x = encrypt(123) // cypher
y = encrypt(456) // cypher
plaintext = 789    // plain
a = x + plaintext // fast
b = x + y // slow
```

## Finite Field Operation

### Modular Reduction

Modular reduction or “mod” is one of the most used operation since in finite field math we basically have to mod everything by some number to make it secure, or mostly $p$. Making this operation fast have tremendous effect on final performance since it’s used basically everywhere.

We tested a few implementation as follow

### Fast Mod

We want to calculate $x \mod n$ where $x < 2n$. In this case, we can compare $x$ with $n$. If $x > n$, subtract $x$ by $n$.

```rust
pub fn modulo_fast<const NB: usize, P: Numeral>(
    x: &RadixCiphertext,
    b: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let len = x.blocks().len();
    let mut x = x.clone();
    let is_gt = server_key.scalar_ge_parallelized(&x, b);
    let to_sub = selector_zero_constant::<NB, _>(b, &is_gt, server_key);
    server_key.sub_assign_parallelized(&mut x, &to_sub);
    server_key.trim_radix_blocks_msb_assign(&mut x, len - NB);
    x
}
```

### General Case Mod

In the `tfhe-rs 0.3` library, there is an implementation of the `div_rem` algorithm that provides a straightforward solution for finding remainders in all cases. However, it is important to note that this implementation is considered to be naive and is known to be one of the slowest operations in the TFHE library.

We can use a algorithm specifically to only reduce $x$ to $x \mod n$, one of that algorithm which is the most simple one is called Barrett reduction. 

### Barrett Reduction

We can simply calculate $x\mod n=x-\lfloor x/n\rfloor \cdot n$ which we can approximate $1/n$ with $m/2^k$, observe that $1/2^k$ is simply right shift by $k$. So $x \mod n$ in this case requires 1 $b\cdot 2$ bits constant mul, 1 $k$-constant bitshift, 1 $b$ bits constant mul, 1 sub, and 1 fast mod.

```rust
// pseudo code
func reduce(a uint) uint {
    q := (a * m) >> k
    a -= q * n
    if a >= n {
        a -= n
    }
    return a
}
```

```rust
// actual implementation
pub fn mod<const NB: usize, P: Numeral>(
    x: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let k = 4 * NB;
    let m_bigint = BigInt::from(2).pow(k as u32) / to_bigint(p);
    let block_to_add = (m_bigint.bits() - 2 * NB as u64 + 1) / 2;
    let m = from_bigint::<U512>(&m_bigint);
    let mut x =
        server_key.extend_radix_with_trivial_zero_blocks_msb(x, NB + block_to_add as usize);
    let mut q = server_key.scalar_mul_parallelized(&x, m);
    server_key.scalar_right_shift_assign_parallelized(&mut q, k as u64);
    server_key.sub_assign_parallelized(&mut x, &server_key.scalar_mul_parallelized(&q, p));
    let len = x.blocks().len();
    server_key.trim_radix_blocks_msb_assign(&mut x, len - (NB + 1));

    modulo_fast::<NB, _>(&x, p, server_key)
}
```

### Mersenne Mod

We saw that Secp256k1 have a prime modulo that match a form of Pseudo-Mersenne prime or prime in the form of $p=2^m - k$ where $k < 2^m/2$ which allow for faster algorithm. One of the algorithm is called Mersenne Reduction that work on $x \mod m$ where $x < n^2$.

We can find $a,b$ from $x= a\cdot 2^m+b$ and then calculate $x \mod n = a\cdot k+b$ and repeat the algorithm until $x < n$. We observe that the algorithm 2 times and then fast mod to get $x \mod n$ without force exit condition.

In this algorithm, we ran 2 pass of Mersenne reduction which ended up in x2 2 $m$-bitshift, x2 1 $b/2$ bits constant mul, and finally 1 fast mod.

```rust
pub fn mod_mersenne<const NB: usize, P: Numeral>(
    x: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let (n, c) = mersenne_coeff_p(p);
    let ceilc = bigint_ilog2_ceil(&c);
    let c_blocks = (c.bits() as usize + 1) / 2;
    let x = server_key.extend_radix_with_trivial_zero_blocks_msb(x, (NB * 2) - x.blocks().len());

    // first pass NB*2 blocks
    let x_mod_p = (|x: &RadixCiphertext| {
        let mut a = server_key.scalar_right_shift_parallelized(x, n as u64);
        let mut b = server_key
            .sub_parallelized(x, &server_key.scalar_left_shift_parallelized(&a, n as u64));

        let len = x.blocks().len();
        // a will be multiplied by c, so it must be at least NB + c_blocks long
        server_key.trim_radix_blocks_msb_assign(&mut a, len - (NB + c_blocks));
        // b must be at least NB long
        server_key.trim_radix_blocks_msb_assign(&mut b, len - NB);
        let ca = server_key.scalar_mul_parallelized(&a, bigint_to_u128(&c));
        server_key.add_parallelized(&ca, &b)
    })(&x);

    // second pass % NB + c_blocks blocks
    let x_mod_p2 = (|x: &RadixCiphertext| {
        let mut a = server_key.scalar_right_shift_parallelized(x, n as u64);
        let mut b = server_key
            .sub_parallelized(x, &server_key.scalar_left_shift_parallelized(&a, n as u64));

        let len = x.blocks().len();
        // a will be multiplied by c, so it must be at least NB + 1 long
        server_key.trim_radix_blocks_msb_assign(&mut a, len - (NB + 1));
        // b must be at least NB long
        server_key.trim_radix_blocks_msb_assign(&mut b, len - NB);
        let ca = server_key.scalar_mul_parallelized(&a, bigint_to_u128(&c));
        server_key.add_parallelized(&ca, &b)
    })(&x_mod_p);

    modulo_fast::<NB, _>(&x_mod_p2, p, server_key)
}
```

### Modular Addition & Modular Multiplication

Our `add_mod` is based on add and fast mod.

```rust
pub fn add_mod<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let mut a_expanded = server_key.extend_radix_with_trivial_zero_blocks_msb(a, 1);
    server_key.add_assign_parallelized(&mut a_expanded, b);
    modulo_fast::<NB, _>(&a_expanded, p, server_key)
}
```

Our `mul_mod` is based on multiplication and Mersenne Barrett Reduction.

```rust
pub fn mul_mod_mersenne<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let mut a_expanded = server_key.extend_radix_with_trivial_zero_blocks_msb(a, NB);
    server_key.mul_assign_parallelized(&mut a_expanded, b);
    mod_mersenne::<NB, _>(&a_expanded, p, server_key)
}
```

### Modular Inverse

Modular inverse $x^{-1}$ is a multiplicative inverse of some value $x$ on $\mod n$ which $x\cdot x^{-1}=1\mod n$ and is considered the slowest function out of all primitive function. But while rarely used, this operation contributed to significant % of final run time. There’re a few way to get Inverse namely

### Fermat's little theorem

We can calculate inverse using $x^{-1} = x^{p-1} \mod p$. While this is a simple algorithm but it require 1 modular exponentiation or $b$ bits double mod, $b$ bits add mod, and $b$ bits mul mod which is extremely expensive.

### Trimmed Extended GCD

This [algorithm](https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm) is faster mainly because it operate over normal number instead of over $\mod p$ and with guarantee that input size will reduce be 1 bit every step, we can trim cypher text to speed up even more. This trimming/bit reducing trick result in over 40% speed up.

```rust
// pseudo code
function extended_gcd(a, b)
    (old_r, r) := (a, b)
    (old_s, s) := (1, 0)
    (old_t, t) := (0, 1)
    
    while r ≠ 0 do
				// each loop will reduce r by at least half (or 1 bit)
        quotient := old_r div r
        (old_r, r) := (r, old_r − quotient × r)
        (old_s, s) := (s, old_s − quotient × s)
        (old_t, t) := (t, old_t − quotient × t)
   
    output "Bézout coefficients:", (old_s, old_t)
    output "greatest common divisor:", old_r
    output "quotients by the gcd:", (t, s)
```

To modify this algorithm to works in homomorphic function. We need to do the $r$ calculation loop without force exit condition, which is $b+1$ round. Also, we need to select correct $r$ when exit condition happened. And finally, the result from this algorithm might be **negative**, so we need to convert it back to **positive**.

```rust
pub fn inverse_mod_trim<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let padded_nb = NB + 1;
    // implement extended euclidean algorithm with trim bit
    // assume a < p. (no check)
    let a = server_key.extend_radix_with_trivial_zero_blocks_msb(&a.clone(), 1);
    let mut r0: RadixCiphertext = server_key.create_trivial_radix(p, padded_nb);
    let mut r1 = a;
    let mut was_done = server_key.create_trivial_radix(0, 1);
    let mut t0 = server_key.create_trivial_radix(0, padded_nb);
    let mut t1: RadixCiphertext = server_key.create_trivial_radix(1, padded_nb);
    let mut inv = server_key.create_trivial_radix(0, padded_nb);
    let mut trim = 0;
    // euclidean algorithm
    // NB/2 best case and NB worst case
    let loop_end = <P as Numeric>::BITS + 1;
    for i in 0..loop_end {
        let _tmr = timer!(Level::Trace; "Inverse Mod", "Bit {}", i);
        // q, r = r0 / r1
        let (mut q, r) = server_key.div_rem_parallelized(&r0, &r1);
        server_key.extend_radix_with_trivial_zero_blocks_msb_assign(&mut q, trim);
        let full_r = server_key.extend_radix_with_trivial_zero_blocks_msb(&r, trim);
        let tmp = t1.clone();
        let qt1 = server_key.mul_parallelized(&t1, &q);
        // t1 = t0 - q * t1
        t1 = server_key.sub_parallelized(&t0, &qt1);
        t0 = tmp;
        // is_done = r =? 0
        // never_done = 1 - is_done
        // was_done = was_done | is_done
        // done_now = is_done & never_done
        let mut done = server_key.scalar_eq_parallelized(&full_r, 0);
        let len = done.blocks().len();
        server_key.trim_radix_blocks_msb_assign(&mut done, len - 1);
        let never_done =
            server_key.sub_parallelized(&server_key.create_trivial_radix(1, 1), &was_done);
        let done_now = server_key.bitand_parallelized(&done, &never_done);
        server_key.bitor_assign_parallelized(&mut was_done, &done);

        let update = selector_zero(&t0, &done_now, server_key);
        server_key.add_assign_parallelized(&mut inv, &update);

        // update values
        if (i % 2 == 0) & (i != 0) {
            r0 = server_key.trim_radix_blocks_msb(&r1.clone(), 1);
            r1 = server_key.trim_radix_blocks_msb(&r.clone(), 1);
            trim += 1;
        } else {
            r0 = r1.clone();
            r1 = r.clone();
        }
    }

    // final result mod p
    // inverse can be **negative**. so we need to add p to make it positive
    server_key.scalar_add_assign_parallelized(&mut inv, p);
    let mut is_gt = server_key.scalar_ge_parallelized(&inv, p);
    server_key.trim_radix_blocks_msb_assign(&mut is_gt, padded_nb - 1);
    let to_sub =
        server_key.mul_parallelized(&server_key.create_trivial_radix(p, padded_nb), &is_gt);
    server_key.sub_assign_parallelized(&mut inv, &to_sub);
    server_key.full_propagate_parallelized(&mut inv);

    inv
}
```

### Binary GCD

There’re also another candidate that we explore but due to it originally have 3 branches, and fast operation (all native operation can be expressed bitshift and constant homomorphic operation).

But finally this algorithm cannot do guaranteed bit trimming so this is slower than trimmed extended GCD.

```python
# python code
def binary_gcd(x, p):
	a = x
	b = p
	u = 1
	v = 0
	i = 0
	while i < bits_length:
	    i += 1
	    s = int(a % 2 == 1)
	    if a % 2 == 1 and a < b:
	        a, b = b, a  # swap
	        u, v = v, u  # swap
	    a = a - b * s
	    u = u - v * s
	    a = int(a / 2)  # sub, shift
	    u = (u * 126) % 251  # mul, mod
	
	print(i)
```

## Elliptic Curve Operation

For elliptic curve operation, we considered using **Affine** and **Projective** **Jacobian** coordinates for group operation calculation.

For illustration, affine coordinates are normal $(x,y)$ pair, but projective coordinates add $z$ component for easier and faster calculation resulted in $(x,y,z)$ triples. Jacobian refers to a type of projective coordinates calculation method. There is more projective coordinates type, for example, **Homogenous** coordinates.

Jacobian triples $(x,y,z)$ represents the affine points $(x/z^2,y/z^3)$ and affine point $(x,y)$ represent the Jacobian triples $(x, y, 1)$

### Group Add

[Adding points](https://en.wikipedia.org/wiki/Elliptic_curve_point_multiplication) in elliptic curve is define as drawing a line through the two points, finding the third intersection point on the curve, and then reflecting it on x axis. Solving for the point can be done in many way with different performance and requirement.

### Affine

Normally in affine coordinates, we can calculate point addition of $(x_1,y_1)$ and $(x_2,y_2)$ by calculating gradient of both point

$$
\lambda =\frac{(y_1−y_2)}{(x_1−x_2)}
$$

Then calculate the final result point directly

$$
x_3=\lambda ^2−x_1−x_2,y_3=\lambda \cdot (x_1−x_2)−y_1
$$

This method seems easy, but it involve 1 modular inverse which is super slow. So we came up with the next trick.

### Jacobian With Jacobian

[Jacobian calculation](https://en.wikibooks.org/wiki/Cryptography/Prime_Curve/Jacobian_Coordinates) benefits from no modular inverse but the points will stay in Projective Jacobian until we convert it back to Affine. This method involves 12 $b$ bits mul mod and 4 $b$ bits add mod.

```python
def jac(X1, Y1, Z1, X2, Y2, Z2):
	U1 = X1*Z2^2
	U2 = X2*Z1^2
	S1 = Y1*Z2^3
	S2 = Y2*Z1^3
	if (U1 == U2)
	 if (S1 != S2)
	   return POINT_AT_INFINITY
	 else
	   return POINT_DOUBLE(X1, Y1, Z1)
	H = U2 - U1
	R = S2 - S1
	X3 = R^2 - H^3 - 2*U1*H^2
	Y3 = R*(U1*H^2 - X3) - S1*H^3
	Z3 = H*Z1*Z2
	return (X3, Y3, Z3)
```

### Jacobian With Affine

By using Jacobian with Jacobian addition and replace $z_2$ with 1, we can reduce the operation used to 8 $b$ bits mul mod and 3 $b$ bits add mod.

We use this algorithm in our scalar multiplication later, In the version of homomorphic form, since there’s no force exit condition, we also need to handle selector if one or both points are zero too.

```rust
pub fn group_projective_add_affine<const NB: usize, P: Numeral>(
    x: &RadixCiphertext,
    y: &RadixCiphertext,
    z: &RadixCiphertext,
    other_x: &RadixCiphertext,
    other_y: &RadixCiphertext,
    other_flag_bit: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext, RadixCiphertext) {
    // z1z1 = z1^2
    let z1z1 = square_mod::<NB, _>(z, p, server_key);
    // u2 = x2*z1z1
    // s2 = y2*z1*z1*z1
    let (u2, s2) = rayon::join(
        || mul_mod::<NB, _>(other_x, &z1z1, p, server_key),
        || {
            mul_mod::<NB, _>(
                other_y,
                &mul_mod::<NB, _>(&z1z1, z, p, server_key),
                p,
                server_key,
            )
        },
    );
    // h = u2 - x1
    let h = sub_mod::<NB, _>(&u2, x, p, server_key);
    // hh = h^2
    let hh = square_mod::<NB, _>(&h, p, server_key);
    // i = 4*hh
    let i = double_mod::<NB, _>(&double_mod::<NB, _>(&hh, p, server_key), p, server_key);
    // j = h*i
    // v = x1*i
    let (j, v) = rayon::join(
        || mul_mod::<NB, _>(&h, &i, p, server_key),
        || mul_mod::<NB, _>(x, &i, p, server_key),
    );
    // r = 2*(s2 - y1)
    let r = double_mod::<NB, _>(&sub_mod::<NB, _>(&s2, y, p, server_key), p, server_key);
    // x3 = r^2 - j - 2*v
    // y3 = r*(v - x3) - 2*y1*j
    // z3 = 2*z1*h
    let ((x3, z3), yj2) = rayon::join(
        || {
            rayon::join(
                || {
                    sub_mod::<NB, _>(
                        &sub_mod::<NB, _>(
                            &square_mod::<NB, _>(&r, p, server_key),
                            &j,
                            p,
                            server_key,
                        ),
                        &double_mod::<NB, _>(&v, p, server_key),
                        p,
                        server_key,
                    )
                },
                || double_mod::<NB, _>(&mul_mod::<NB, _>(z, &h, p, server_key), p, server_key),
            )
        },
        || mul_mod::<NB, _>(y, &double_mod::<NB, _>(&j, p, server_key), p, server_key),
    );
    let y3 = sub_mod::<NB, _>(
        &mul_mod::<NB, _>(&r, &sub_mod::<NB, _>(&v, &x3, p, server_key), p, server_key),
        &yj2,
        p,
        server_key,
    );

    // z1'/z0' 0  1
    //    0    x' x1
    //    1    x0 x0
    // x'' =  x' * is_z0_z1_non_zero + (x0 + x1) * not_is_z0_z1_non_zero
    // y'' =  y' * is_z0_z1_non_zero + (y0 + y1) * not_is_z0_z1_non_zero
    // z'' =  z' * is_z0_z1_non_zero + (z0 + z1) * not_is_z0_z1_non_zero
    let (mut is_z0_non_zero, mut is_z1_non_zero) = rayon::join(
        || server_key.scalar_ne_parallelized(z, 0),
        || server_key.scalar_ne_parallelized(other_flag_bit, 0),
    );
    server_key.trim_radix_blocks_msb_assign(&mut is_z0_non_zero, NB - 1);
    server_key.trim_radix_blocks_msb_assign(&mut is_z1_non_zero, NB - 1);
    let is_z0_z1_non_zero = server_key.bitand_parallelized(&is_z0_non_zero, &is_z1_non_zero);
    let not_is_z0_z1_non_zero =
        server_key.sub_parallelized(&server_key.create_trivial_radix(1, 1), &is_z0_z1_non_zero);

    let (((xp1, xp2), (yp1, yp2)), (zp1, zp2)) = rayon::join(
        || {
            rayon::join(
                || {
                    rayon::join(
                        || server_key.mul_parallelized(&x3, &is_z0_z1_non_zero),
                        || {
                            server_key.mul_parallelized(
                                &server_key.add_parallelized(x, other_x),
                                &not_is_z0_z1_non_zero,
                            )
                        },
                    )
                },
                || {
                    rayon::join(
                        || server_key.mul_parallelized(&y3, &is_z0_z1_non_zero),
                        || {
                            server_key.mul_parallelized(
                                &server_key.add_parallelized(y, other_y),
                                &not_is_z0_z1_non_zero,
                            )
                        },
                    )
                },
            )
        },
        || {
            rayon::join(
                || server_key.mul_parallelized(&z3, &is_z0_z1_non_zero),
                || {
                    server_key.mul_parallelized(
                        &server_key.add_parallelized(z, other_flag_bit),
                        &not_is_z0_z1_non_zero,
                    )
                },
            )
        },
    );

    let ((x_prime, y_prime), z_prime) = rayon::join(
        || {
            rayon::join(
                || server_key.add_parallelized(&xp1, &xp2),
                || server_key.add_parallelized(&yp1, &yp2),
            )
        },
        || server_key.add_parallelized(&zp1, &zp2),
    );

    (x_prime, y_prime, z_prime)
}
```

### Group Double

Adding points but adding with itself is a special case where tangent line is used instead of a standard secant line. This algorithm is not needed to be homomorphic in this work.

### Group Scalar Multiplication

This operation only used once but contribute to almost half of the final run time. Scalar multiplication $xG$ basically adding $G$ to itself $x$ times. The main native algorithm that used to operate this efficiently is double-add method

### Double-And-Add

Basically this algorithm iterate through each bits of $x$, double the $G$ in between. If the bit is 1, add that current $G$ to the result.

```python
# python code
def double_and_add(G, x):
	result = identity_point
	for b in bits(x):
	    if b == 1:
	        result = result + G // add if current bit = 1
			P = double(G)
	return result
```

In the homomorphic variant, we calculate the bits by bitand the current bits by 1, and update the current bits by right bitshift by 1 every iteration. So in total it requires $b$ $b$-bits bitand, bitshift, projective addition, projective double, and $3b$ $b$-bits multiplication.

```rust
pub fn group_projective_scalar_mul<const NB: usize, P: Numeral>(
    x: &RadixCiphertext,
    y: &RadixCiphertext,
    z: &RadixCiphertext,
    scalar: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext, RadixCiphertext) {
    let mut tmp_x = x.clone();
    let mut tmp_y = y.clone();
    let mut tmp_z = z.clone();
    let mut scalar = scalar.clone();
    let mut res_x = server_key.create_trivial_radix(0, NB);
    let mut res_y = server_key.create_trivial_radix(0, NB);
    let mut res_z = server_key.create_trivial_radix(0, NB);

    for _i in 0..<P as Numeric>::BITS {
        let (mut bit, new_scalar) = rayon::join(
            || server_key.scalar_bitand_parallelized(&scalar, 1),
            || server_key.scalar_right_shift_parallelized(&scalar, 1),
        );
        server_key.trim_radix_blocks_msb_assign(&mut bit, NB - 1);
        scalar = new_scalar;
        ((res_x, res_y, res_z), (tmp_x, tmp_y, tmp_z)) = rayon::join(
            || {
                let ((x_to_add, y_to_add), z_to_add) = rayon::join(
                    || {
                        rayon::join(
                            || server_key.mul_parallelized(&tmp_x, &bit),
                            || server_key.mul_parallelized(&tmp_y, &bit),
                        )
                    },
                    || server_key.mul_parallelized(&tmp_z, &bit),
                );
                group_projective_add_projective::<NB, _>(
                    &res_x, &res_y, &res_z, &x_to_add, &y_to_add, &z_to_add, p, server_key,
                )
            },
            || group_projective_double::<NB, _>(&tmp_x, &tmp_y, &tmp_z, p, server_key),
        );
    }

    (res_x, res_y, res_z)
}
```

### Sliding Window

In our case, we need scalar multiplication $kG$ for scalar $k$ on generator point $G$. which mean that we can compute all any $xG$ efficiently. and then select result based on scalar input $k$.

However, the amount of possible outcome scale exponentially by the amount of scalar total bits. Which mean that it not practical for 256 bits which need $2^{256}$ outcome, but we can do [sliding window](https://www.notion.so/tutorial-ffdd7366b07441368b8e5b966757c0b2?pvs=21) on the value to select some of the bits from the scalar and make the selector less computation intensive.

So basically, we select first $W$ bits from the scalar, precompute $2^{W-1}$ points by doubling $G$ in native function, make $2^{W-1}$ selector and select 1 point from those precomputed points, add the selected point to the result, traverse further by $W$ bits, then repeat until bits are exhausted.

```rust
// pseudo code
Q ← 0
for i from m downto 0 do
    if di = 0 then
        Q ← point_double(Q)
    else 
        t ← extract j (up to w − 1) additional bits from d (including di)
        i ← i − j
        if j < w then
            Perform double-and-add using t 
            return Q
        else 
            Q ← point_double_repeat(Q, w)
            Q ← point_add(Q, tP)
return Q
```

In the homomorphic variants, we can reuse the bits selection from the double-and-add method, calculate not-bit, bitand each bit together to select the value. In the selection process we can use tree-based calculation method to parallelize the multiplication and addition process.

```rust
pub fn group_projective_scalar_mul_constant_windowed<
    const W: usize,
    const NB: usize,
    P: Numeral,
>(
    x: P,
    y: P,
    scalar: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext, RadixCiphertext) {
    let mut tmp_x = x;
    let mut tmp_y = y;
    let mut scalar = scalar.clone();
    let mut res_x = server_key.create_trivial_radix(0, NB);
    let mut res_y = server_key.create_trivial_radix(0, NB);
    let mut res_z = server_key.create_trivial_radix(0, NB);

    // take W bits at a time
    // for each bit, we have a precomputed points of 2^W - 1 points
    // take the bit, and use it to select the point
    // add the point to the result
    // then double the temporary point W times
    let mut i = 0;
    while i < <P as Numeric>::BITS {
        let chunk_size = match i + W > <P as Numeric>::BITS {
            true => <P as Numeric>::BITS - i,
            false => W,
        };
        let _ic = i..i + chunk_size;
        i += chunk_size;

        let _tmr = stimer!(Level::Info; "Scalar Mul", "Bits {:?}", _ic);
        let cal_bits_tmr = timer!(Level::Debug; "Calculating bits");
        // get the next W bits
        let mut tmp_bits = vec![
            (
                server_key.create_trivial_radix(0, NB),
                server_key.create_trivial_radix(0, NB),
            );
            chunk_size
        ];
        (0..chunk_size)
            .into_par_iter()
            .map(|i| {
                let shifted = server_key.scalar_right_shift_parallelized(&scalar, i as u64);
                let mut bit = server_key.scalar_bitand_parallelized(&shifted, 1);
                server_key.trim_radix_blocks_msb_assign(&mut bit, NB - 1);
                (
                    server_key.sub_parallelized(&server_key.create_trivial_radix(P::ONE, 1), &bit),
                    bit,
                )
            })
            .collect_into_vec(&mut tmp_bits);
        let mut bits = vec![];
        let mut not_bits = vec![];
        for (not_bit, bit) in tmp_bits {
            not_bits.push(not_bit);
            bits.push(bit);
        }
        server_key.scalar_right_shift_assign_parallelized(&mut scalar, chunk_size as u64);
        drop(cal_bits_tmr);

        // get the precomputed values
        let mut points = vec![(P::ZERO, P::ZERO)];
        let tmp = (tmp_x, tmp_y);
        for _ in 1..2usize.pow(chunk_size as u32) {
            points.push((tmp_x, tmp_y));
            // points are stored in tmp
            (tmp_x, tmp_y) = {
                let (tmp_x_new, temp_y_new, temp_z_new) =
                    group_projective_add_affine_native(tmp_x, tmp_y, P::ONE, tmp.0, tmp.1, p);
                group_projective_into_affine_native(tmp_x_new, temp_y_new, temp_z_new, p)
            };
        }

        // select the points
        let sel_tmr = timer!(Level::Debug; "Selecting points", "Points {}", points.len() - 1);
        let mut points_to_add = vec![
            (
                server_key.create_trivial_radix(0, NB),
                server_key.create_trivial_radix(0, NB)
            );
            2usize.pow(chunk_size as u32) - 1
        ];
        points
            .into_par_iter()
            .enumerate()
            .take(2usize.pow(chunk_size as u32))
            .skip(1)
            .map(|(i, point)| {
                let bits = (0..chunk_size)
                    .map(|j| match i & 2usize.pow(j as u32) == 0 {
                        true => not_bits[j].clone(),
                        false => bits[j].clone(),
                    })
                    .collect::<Vec<_>>();
                let selected_bit =
                    parallel_fn(&bits, |b0, b1| server_key.bitand_parallelized(b0, b1));
                rayon::join(
                    || selector_zero_constant::<NB, _>(point.0, &selected_bit, server_key),
                    || selector_zero_constant::<NB, _>(point.1, &selected_bit, server_key),
                )
            })
            .collect_into_vec(&mut points_to_add);
        let selected_point = parallel_fn(&points_to_add, |p0, p1| {
            rayon::join(
                || server_key.add_parallelized(&p0.0, &p1.0),
                || server_key.add_parallelized(&p0.1, &p1.1),
            )
        });
        drop(sel_tmr);

        // check if all bits are not zero for flag bit
        let kary_or_tmr = timer!(Level::Debug; "Kary or");
        let all_not_zero = parallel_fn(&bits, |b0, b1| server_key.bitor_parallelized(b0, b1));
        drop(kary_or_tmr);

        // add the point
        (res_x, res_y, res_z) = group_projective_add_affine::<NB, _>(
            &res_x,
            &res_y,
            &res_z,
            &selected_point.0,
            &selected_point.1,
            &all_not_zero,
            p,
            server_key,
        );
    }

    (res_x, res_y, res_z)
}
```

We can estimate the total cost of different window size with this equation

$$
⁍ 
$$

Where

- $w$ is the window size, e.g 4, 5, 6
- $T_S$ is the computation time per selector
- $T_K$ is the group projective add affine time
- $B$ is the total bits

We observed that based on 64 cores machine, the optimal window size are 6 or 7.

![window.png](res/window.png)

### Group Transformation

Simply convert projective Jacobian form back to affine form. This requires 1 inverse mod and 2 mul mod.

```rust
pub fn group_projective_into_affine_inv<const NB: usize, P: Numeral>(
    x: &RadixCiphertext,
    y: &RadixCiphertext,
    z_inv: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext) {
    let z_inv2 = square_mod::<NB, _>(z_inv, p, server_key);
    let z_inv3 = mul_mod::<NB, _>(&z_inv2, z_inv, p, server_key);

    rayon::join(
        || mul_mod::<NB, _>(x, &z_inv2, p, server_key),
        || mul_mod::<NB, _>(y, &z_inv3, p, server_key),
    )
}
```

## ECDSA

Finally, we can put everything together. in our implementation we assume that private key $sk$, message, and nonce $k$ are provided.

We can calculate $z^{-1}$ and $k^{-1}$ in parallel to save sometime. 

```rust
pub fn ecdsa_sign<const NB: usize, P: Numeral>(
    sk: &RadixCiphertext,
    k: &RadixCiphertext,
    message: P,
    generator: (P, P),
    q_modulo: P,
    r_modulo: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext) {
    // (x, y) = k * G
    let (x_proj, y_proj, z_proj) = group_projective_scalar_mul_constant_windowed::<WINDOW, NB, _>(
        generator.0,
        generator.1,
        k,
        q_modulo,
        server_key,
    );
    let (z_inv, k_inv) = rayon::join(
        || inverse_mod::<NB, _>(&z_proj, q_modulo, server_key),
        || inverse_mod::<NB, _>(k, r_modulo, server_key),
    );
    let (x, y) =
        group_projective_into_affine_inv::<NB, _>(&x_proj, &y_proj, &z_inv, q_modulo, server_key);
    // r = x
    // s = k^-1 * (m + r * sk)
    let r = if q_modulo > r_modulo && q_modulo <= P::TWO * r_modulo {
        modulo_fast::<NB, _>(&x, r_modulo, server_key)
    } else {
        mod_mersenne::<NB, _>(&x, r_modulo, server_key)
    };
    let mrsk = add_mod::<NB, _>(
        &server_key.create_trivial_radix(message, NB),
        &mul_mod::<NB, _>(&r, sk, r_modulo, server_key),
        r_modulo,
        server_key,
    );
    let s = mul_mod::<NB, _>(&k_inv, &mrsk, r_modulo, server_key);
    (r, s)
}
```

Additionally with the same building blocks, we can also implement the ECDSA verification function.

```rust
pub fn ecdsa_verify<const NB: usize, P: Numeral>(
    public_key: (P, P),
    signature: (RadixCiphertext, RadixCiphertext),
    message: P,
    generator: (P, P),
    q_modulo: P,
    r_modulo: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    // s^-1
    let s_inv = inverse_mod::<NB, _>(&signature.1, r_modulo, server_key);
    // u1 = m * s^-1
    let u1 = mul_mod::<NB, _>(
        &s_inv,
        &server_key.create_trivial_radix(message, NB),
        r_modulo,
        server_key,
    );
    // u2 = r * s^-1
    let u2 = mul_mod::<NB, _>(&s_inv, &signature.0, r_modulo, server_key);
    // (x, y) = u1 * G + u2 * Q
    let (x_proj_1, y_proj_1, z_proj_1) = group_projective_scalar_mul_constant_windowed::<
        WINDOW,
        NB,
        _,
    >(generator.0, generator.1, &u1, q_modulo, server_key);
    let (x_proj_2, y_proj_2, z_proj_2) = group_projective_scalar_mul_constant_windowed::<
        WINDOW,
        NB,
        _,
    >(
        public_key.0, public_key.1, &u2, q_modulo, server_key
    );
    let (x_proj, y_proj, z_proj) = group_projective_add_projective::<NB, _>(
        &x_proj_1, &y_proj_1, &z_proj_1, &x_proj_2, &y_proj_2, &z_proj_2, q_modulo, server_key,
    );
    let is_z_zero = server_key.scalar_eq_parallelized(&z_proj, 0);
    let (x, _y) =
        group_projective_into_affine::<NB, _>(&x_proj, &y_proj, &z_proj, q_modulo, server_key);
    let x_mod_scalar = modulo_fast::<NB, _>(&x, r_modulo, server_key);
    let is_x_eq_r = server_key.eq_parallelized(&x_mod_scalar, &signature.0);

    // valid if z != 0 && x == r
    server_key.bitand_parallelized(&is_z_zero, &is_x_eq_r)
} 
```

## Run this code

### setup

- this code will run for 1-2 day on 64 cores machine. cloud VM and tmux is recommended but not required
- [rust](https://www.rust-lang.org/tools/install)

### run

```bash
git clone https://github.com/Tetration-Lab/TFHECDSA
cd TFHECDSA
cargo run --release --example ecdsa
```

### log config

log level, there’s 3 level `Info`, `Debug`, `Traces` where default is `Debug` to change this use use `RUST_LOG="level"`  

### navigation

```bash
src
├── ecdsa.rs // main entry point here
├── helper.rs 
├── lib.rs
├── main.rs
├── numeral.rs
├── ops 
│   ├── group_jacobian.rs
│   ├── mersenne.rs
│   ├── native.rs
│   ├── primitive.rs
│   └── secp256k1.rs
├── ops.rs
└── stats.rs
examples
└── ecdsa.rs // see how to use here
```