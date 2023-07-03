use std::time::Instant;

use tfhe::{
    core_crypto::prelude::{Numeric, UnsignedInteger},
    integer::ClientKey,
    integer::{
        block_decomposition::{DecomposableInto, RecomposableFrom},
        IntegerCiphertext, RadixCiphertext, ServerKey,
    },
};

// a^-1 mod p where a*a^-1 = 1 mod p
pub fn inverse_mod<
    const NB: usize,
    P: DecomposableInto<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    a: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    // implement extended euclidean algorithm
    // assume a < p. (no check)
    let a = a.clone();
    let mut r0 = server_key.create_trivial_radix(p, NB);
    let mut r1 = a.clone();
    let mut was_done = server_key.create_trivial_radix(0, 1);
    let mut t0 = server_key.create_trivial_radix(0, NB);
    let mut t1 = server_key.create_trivial_radix(1, NB);
    let mut inv = server_key.create_trivial_radix(0, NB);

    // euclidean algorithm
    // NB/2 best case and NB worst case
    for _ in 0..<P as Numeric>::BITS / 2 {
        let (q, r) = server_key.smart_div_rem_parallelized(&mut r0.clone(), &mut r1.clone());
        let tmp = t1.clone();
        t1 = server_key.smart_sub_parallelized(
            &mut t0.clone(),
            &mut server_key.smart_mul_parallelized(&mut q.clone(), &mut t1.clone()),
        );
        t0 = tmp;
        // is_done = r =? 0
        // never_done = 1 - is_done
        // was_done = was_done | is_done
        // done_now = is_done & never_done
        let mut done = server_key.smart_scalar_eq_parallelized(&mut r.clone(), 0);
        let mut never_done = server_key
            .smart_sub_parallelized(&mut server_key.create_trivial_radix(1, 1), &mut was_done);
        let mut done_now = server_key.smart_bitand_parallelized(&mut done, &mut never_done);
        server_key.smart_bitor_assign_parallelized(&mut was_done, &mut done);
        let mut update = server_key.smart_mul_parallelized(&mut done_now, &mut t0);
        server_key.smart_add_assign_parallelized(&mut inv, &mut update);
        // update values
        r0 = r1;
        r1 = r;
    }

    // final result
    let mut tmp = server_key.smart_scalar_add_parallelized(&mut inv, p);
    let mut is_gt = server_key.smart_scalar_ge_parallelized(&mut tmp, p);
    server_key.trim_radix_blocks_msb_assign(&mut is_gt, NB - 1);
    let mut to_sub =
        server_key.smart_mul_parallelized(&mut server_key.create_trivial_radix(p, NB), &mut is_gt);
    server_key.smart_sub_assign_parallelized(&mut tmp, &mut to_sub);
    tmp
}

use crate::helper::format;

/// a + b mod p
pub fn add_mod<const NB: usize, P: DecomposableInto<u64> + DecomposableInto<u8> + Copy + Sync>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    // assume large p and a,b < p
    let mut a_expanded = server_key.extend_radix_with_trivial_zero_blocks_msb(a, 1);
    server_key.smart_add_assign_parallelized(&mut a_expanded, &mut b.clone());
    let mut is_gt = server_key.smart_scalar_gt_parallelized(&mut a_expanded, p);
    server_key.trim_radix_blocks_msb_assign(&mut is_gt, NB - 1);
    let mut to_sub =
        server_key.smart_mul_parallelized(&mut server_key.create_trivial_radix(p, NB), &mut is_gt);
    server_key.smart_sub_assign_parallelized(&mut a_expanded, &mut to_sub);
    server_key.trim_radix_blocks_msb_assign(&mut a_expanded, 1);
    a_expanded
}

/// a - b mod p
pub fn sub_mod<const NB: usize, P: DecomposableInto<u64> + DecomposableInto<u8> + Copy + Sync>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let mut is_gt = server_key.smart_gt_parallelized(&mut b.clone(), &mut a.clone());
    let mut to_add =
        server_key.smart_mul_parallelized(&mut server_key.create_trivial_radix(p, NB), &mut is_gt);
    let mut a_expanded = server_key.extend_radix_with_trivial_zero_blocks_msb(a, 1);
    server_key.trim_radix_blocks_msb_assign(&mut is_gt, NB - 1);
    server_key.smart_add_assign_parallelized(&mut a_expanded, &mut to_add);
    server_key.smart_sub_assign_parallelized(&mut a_expanded, &mut b.clone());
    server_key.trim_radix_blocks_msb_assign(&mut a_expanded, 1);
    a_expanded
}

/// a * b mod p
pub fn mul_mod_bitwise<
    const NB: usize,
    P: DecomposableInto<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    // assume large p and a,b < p
    let now = Instant::now();

    let mut res = server_key.create_trivial_radix(0u64, NB);
    let mut a_tmp = a.clone();
    let b_tmp = b.clone();
    let (mut b_next_tmp, mut bit) = rayon::join(
        || server_key.scalar_right_shift_parallelized(&b_tmp, 1),
        || server_key.smart_scalar_bitand_parallelized(&mut b_tmp.clone(), 1),
    );
    server_key.trim_radix_blocks_msb_assign(&mut bit, NB - 1);
    #[allow(clippy::redundant_clone)]
    let mut to_add_later = res.clone();

    println!("initial cost - {}ms", now.elapsed().as_millis());

    for i in 0..<P as Numeric>::BITS {
        let now = Instant::now();

        ((b_next_tmp, a_tmp), (bit, (res, to_add_later))) = rayon::join(
            || {
                rayon::join(
                    || server_key.scalar_right_shift_parallelized(&b_next_tmp, 1),
                    || add_mod::<NB, _>(&a_tmp, &a_tmp, p, server_key),
                )
            },
            || {
                rayon::join(
                    || {
                        let mut bit =
                            server_key.smart_scalar_bitand_parallelized(&mut b_next_tmp.clone(), 1);
                        server_key.trim_radix_blocks_msb_assign(&mut bit, NB - 1);
                        bit
                    },
                    || {
                        rayon::join(
                            || add_mod::<NB, _>(&res, &to_add_later, p, server_key),
                            || server_key.smart_mul_parallelized(&mut a_tmp.clone(), &mut bit),
                        )
                    },
                )
            },
        );

        println!("time used for bit {i} - {}s", now.elapsed().as_secs());
    }

    add_mod::<NB, _>(&res, &to_add_later, p, server_key)
}

/// a * b mod p
pub fn mul_mod<const NB: usize, P: DecomposableInto<u64> + DecomposableInto<u8> + Copy + Sync>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    // assume large p and a,b < p
    let mut a_expanded = server_key.extend_radix_with_trivial_zero_blocks_msb(a, NB);
    server_key.smart_mul_assign_parallelized(&mut a_expanded, &mut b.clone());
    let (_q, mut r) = server_key.smart_div_rem_parallelized(
        &mut a_expanded,
        &mut server_key.create_trivial_radix(p, NB * 2),
    );
    server_key.trim_radix_blocks_msb_assign(&mut r, NB);
    r
}

/// a * b mod p where b is a constant
/// slower than 12 `add_mod`
pub fn mul_mod_constant<
    const NB: usize,
    P: DecomposableInto<u64> + DecomposableInto<u8> + Copy + Sync,
    S: DecomposableInto<u8> + UnsignedInteger,
>(
    a: &RadixCiphertext,
    b: S,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    // assume large p and a,b < p
    let mut a_expanded = server_key.extend_radix_with_trivial_zero_blocks_msb(a, NB);
    server_key.smart_scalar_mul_assign_parallelized(&mut a_expanded, b);
    let (_q, mut r) = server_key.smart_div_rem_parallelized(
        &mut a_expanded,
        &mut server_key.create_trivial_radix(p, NB * 2),
    );
    server_key.trim_radix_blocks_msb_assign(&mut r, NB);
    r
}

/// a^2 mod p
pub fn square_mod<
    const NB: usize,
    P: DecomposableInto<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    a: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    mul_mod::<NB, _>(a, a, p, server_key)
}

/// a*2 mod p
pub fn double_mod<
    const NB: usize,
    P: DecomposableInto<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    a: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    add_mod::<NB, _>(a, a, p, server_key)
}

pub fn group_projective_double<
    const NB: usize,
    P: DecomposableInto<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    x: &RadixCiphertext,
    y: &RadixCiphertext,
    z: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext, RadixCiphertext) {
    let timer = Instant::now();
    // y^2 = x^3+ax+b
    // in case a == 0

    // a = x^2
    // b = y^2
    let (a, b) = rayon::join(
        || square_mod::<NB, _>(x, p, server_key),
        || square_mod::<NB, _>(y, p, server_key),
    );
    println!("a,b - {:.2}s", timer.elapsed().as_secs_f64());
    // c = b^2
    let c = square_mod::<NB, _>(&b, p, server_key);
    println!("c - {:.2}s", timer.elapsed().as_secs_f64());
    // d = 2(2(x + b) - (a + c))
    let (xb2, ac) = rayon::join(
        || double_mod::<NB, _>(&add_mod::<NB, _>(x, &b, p, server_key), p, server_key),
        || add_mod::<NB, _>(&a, &c, p, server_key),
    );
    let d = double_mod::<NB, _>(
        &sub_mod::<NB, _>(
            &xb2, // 2(x + b)
            &ac,  // (a + c)
            p, server_key,
        ),
        p,
        server_key,
    );
    println!("d - {:.2}s", timer.elapsed().as_secs_f64());
    // e = 3a
    let e = add_mod::<NB, _>(&double_mod::<NB, _>(&a, p, server_key), &a, p, server_key);
    println!("e - {:.2}s", timer.elapsed().as_secs_f64());
    // f = e^2
    let f = square_mod::<NB, _>(&e, p, server_key);
    println!("f - {:.2}s", timer.elapsed().as_secs_f64());
    // z' = 2 * y1 * z1
    // x' = f - 2d
    let (z_prime, x_prime) = rayon::join(
        || double_mod::<NB, _>(&mul_mod::<NB, _>(y, z, p, server_key), p, server_key),
        || sub_mod::<NB, _>(&f, &double_mod::<NB, _>(&d, p, server_key), p, server_key),
    );
    println!("z',x' - {:.2}s", timer.elapsed().as_secs_f64());
    // y' = e(d - x') - 8c
    let (edx, c8) = rayon::join(
        || {
            mul_mod::<NB, _>(
                &e,
                &sub_mod::<NB, _>(&d, &x_prime, p, server_key),
                p,
                server_key,
            )
        },
        || {
            let enc_2c = double_mod::<NB, _>(&c, p, server_key);
            let enc_4c = double_mod::<NB, _>(&enc_2c, p, server_key);
            double_mod::<NB, _>(&enc_4c, p, server_key)
        },
    );
    let y_prime = sub_mod::<NB, _>(&edx, &c8, p, server_key);
    println!("y' - {:.2}s", timer.elapsed().as_secs_f64());

    (x_prime, y_prime, z_prime)
}

pub fn group_projective_add_affine<
    const NB: usize,
    P: DecomposableInto<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    x1: &RadixCiphertext,
    y1: &RadixCiphertext,
    z1: &RadixCiphertext,
    x2: &RadixCiphertext,
    y2: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext, RadixCiphertext) {
    // z1z1 = z1^2
    // y2z1 = y2 * z1
    let (z1z1, y2z1) = rayon::join(
        || square_mod::<NB, _>(z1, p, server_key),
        || mul_mod::<NB, _>(y2, z1, p, server_key),
    );
    // u2 = x2 * z1z1
    // s2 = y2z1 * z1z1
    let (u2, s2) = rayon::join(
        || mul_mod::<NB, _>(x2, &z1z1, p, server_key),
        || mul_mod::<NB, _>(&y2z1, &z1z1, p, server_key),
    );
    // h = u2 - x1
    let h = sub_mod::<NB, _>(&u2, x1, p, server_key);
    // hh = h^2
    let hh = square_mod::<NB, _>(&h, p, server_key);
    // i = 4 * hh
    let i = double_mod::<NB, _>(&double_mod::<NB, _>(&hh, p, server_key), p, server_key);
    // j = h * i
    let j = mul_mod::<NB, _>(&h, &i, p, server_key);
    // r = 2 * (s2 - y1)
    let r = double_mod::<NB, _>(&sub_mod::<NB, _>(&s2, y1, p, server_key), p, server_key);
    // v = x1 * i
    let v = mul_mod::<NB, _>(x1, &i, p, server_key);
    // x3 = r^2 - j - 2 * v
    let (r2, j2v) = rayon::join(
        || square_mod::<NB, _>(&r, p, server_key),
        || {
            let j2 = double_mod::<NB, _>(&j, p, server_key);
            sub_mod::<NB, _>(&j2, &double_mod::<NB, _>(&v, p, server_key), p, server_key)
        },
    );
    let x3 = sub_mod::<NB, _>(&r2, &j2v, p, server_key);
    // y3 = r * (v - x3) - 2 * y1 * j
    let (ryx3, y1j2) = rayon::join(
        || mul_mod::<NB, _>(&r, &sub_mod::<NB, _>(&v, &x3, p, server_key), p, server_key),
        || double_mod::<NB, _>(&mul_mod::<NB, _>(y1, &j, p, server_key), p, server_key),
    );
    let y3 = sub_mod::<NB, _>(&ryx3, &y1j2, p, server_key);
    // z3 = (z1 + h)^2 - (z1z1 + hh)
    let (z1h2, z1z1hh) = rayon::join(
        || square_mod::<NB, _>(&add_mod::<NB, _>(z1, &h, p, server_key), p, server_key),
        || add_mod::<NB, _>(z1, &h, p, server_key),
    );
    let z3 = sub_mod::<NB, _>(&z1h2, &z1z1hh, p, server_key);

    (x3, y3, z3)
}

/// a^b mod p
pub fn pow_mod<
    const NB: usize,
    P: DecomposableInto<u64> + RecomposableFrom<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let mut res = server_key.create_trivial_radix(1, NB);
    let mut base = a.clone();
    let mut exponent = b.clone();

    for i in 0..<P as Numeric>::BITS {
        let timer = Instant::now();
        (res, (exponent, base)) = rayon::join(
            || {
                let mut bit = server_key.scalar_bitand_parallelized(&exponent, 1);
                // The line below breaks subtraction
                //server_key.trim_radix_blocks_msb_assign(&mut bit, NB - 1);
                // tmp = bit == 1 ? base : 1;
                // tmp = base * bit + 1 - bit
                let mut tmp = server_key.smart_mul_parallelized(
                    &mut base.clone(),
                    &mut server_key.trim_radix_blocks_msb(&bit, NB - 1),
                );
                server_key.smart_scalar_add_assign_parallelized(&mut tmp, 1);
                server_key.smart_sub_assign_parallelized(&mut tmp, &mut bit);
                mul_mod::<NB, _>(&res, &tmp, p, server_key)
            },
            || {
                rayon::join(
                    || server_key.scalar_right_shift_parallelized(&exponent, 1),
                    || square_mod::<NB, _>(&base, p, server_key),
                )
            },
        );
        println!("pow mod bit {i} took {:.2}s", timer.elapsed().as_secs_f32());
    }

    res
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use tfhe::{
        core_crypto::prelude::Numeric,
        integer::{keycache::IntegerKeyCache, IntegerCiphertext, RadixCiphertext},
        shortint::prelude::PARAM_MESSAGE_2_CARRY_2,
    };

    use crate::{
        helper::format,
        ops::{double_mod, mul_mod_constant},
    };

    #[test]
    #[ignore = "bench"]
    fn bench_8a() {
        let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);

        const NUM_BLOCK: usize = 8;
        type Integer = u16;
        let p: Integer = 13841;
        let a: Integer = 13820;
        let enc_a = client_key.encrypt_radix(a, NUM_BLOCK);

        let timer = Instant::now();
        let enc_2a = double_mod::<NUM_BLOCK, _>(&enc_a, p, &server_key);
        let enc_4a = double_mod::<NUM_BLOCK, _>(&enc_2a, p, &server_key);
        let _enc_8a = double_mod::<NUM_BLOCK, _>(&enc_4a, p, &server_key);
        println!("8a using addition - {}s", timer.elapsed().as_secs());

        let timer = Instant::now();
        let _enc_8a_mul = mul_mod_constant::<NUM_BLOCK, _, _>(&enc_a, 8u8, p, &server_key);
        println!("8a using multiplication - {}s", timer.elapsed().as_secs());
    }

    #[test]
    fn correct_pow_mod() {
        type Integer = u128;
        let p: Integer = 251;
        let x: Integer = 199;
        let y: Integer = 81;

        let mut base = x;
        let mut exponent = y;
        let mut res = <Integer as Numeric>::ONE;

        for _ in 0..<Integer as Numeric>::BITS {
            let bit = exponent & 1;
            let tmp = bit * base + 1 - bit;
            res = res * tmp % p;
            base = (base * base) % p;
            exponent >>= 1;
        }

        println!("Result: {}", res);
    }
}
