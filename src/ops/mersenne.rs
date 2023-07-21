#![allow(clippy::redundant_closure_call)]

use std::time::Instant;

use logging_timer::time;
use num_bigint::BigInt;
use rand::Rng;
use tfhe::integer::{
    block_decomposition::{DecomposableInto, RecomposableFrom},
    IntegerCiphertext, RadixCiphertext, ServerKey, U512,
};

use crate::{
    helper::{bigint_ilog2_ceil, bigint_to_u128, from_bigint, to_bigint},
    numeral::Numeral,
};

use super::{modulo_div_rem, modulo_fast};

/// Calculate n, m, p from coeff
/// `coeff` in the form of p = 2^n_0 - 2^n_1 - ... - 2^n_{k-1} - n_k
/// `c` in the form of c = 2^n_0 - p
/// `c` must be in range 0 <= c <= 2^floor(n/2)
#[inline(always)]
pub fn mersenne_coeff(coeff: &[u32]) -> (u32, BigInt, BigInt, BigInt) {
    assert!(coeff.len() > 1);
    let len = coeff.len();
    let n = coeff[0];
    let p = coeff[1..len - 1]
        .iter()
        .fold(BigInt::from(2).pow(coeff[0]), |acc, b| {
            acc - BigInt::from(2).pow(*b)
        })
        - coeff[len - 1];
    let q = BigInt::from(2).pow(n);
    let c = &q - &p;

    (n, p, q, c)
}

#[inline(always)]
/// Calculate n, c from p
/// `c` must be in range 0 <= c <= 2^floor(n/2)
pub fn mersenne_coeff_p<P: Numeral>(p: P) -> (u32, BigInt) {
    let pb = to_bigint(p);
    let n = bigint_ilog2_ceil(&pb);
    let c = (BigInt::from(1) << n) - &pb;

    (n, c)
}

/// Calculate x mod p^2 mod p
/// `coeff` in the form of p = 2^n_0 - 2^n_1 - ... - 2^n_{k-1} - n_k
pub fn mersenne_mod_native<P: Numeral>(x: P, p: P) -> P {
    let (n, c) = mersenne_coeff_p(p);
    let x = to_bigint(x);
    let p = to_bigint(p);

    // x = a*2^n + b
    let a = &x >> n;
    let b = &x - (&a << n);

    // x % p = (a*2^n + b) % p = c * a + b % p
    let x_mod_p = &a * &c + &b;

    let a = &x_mod_p >> n;
    let b = &x_mod_p - (&a << n);

    // x % p = (a*2^n + b) % p = c * a + b % p
    let x_mod_p = &a * &c + &b;

    from_bigint(&if x_mod_p >= p { x_mod_p - p } else { x_mod_p })
}

/// Calculate x mod p^2 mod p
#[time("trace", "Modulus Reduction Mersenne+Barrett")]
pub fn mod_mersenne<const NB: usize, P: Numeral>(
    x: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let (n, c) = mersenne_coeff_p(p);
    let ceilc = bigint_ilog2_ceil(&c);
    if ceilc >= n / 2 {
        let k = 4 * NB;
        let m_bigint = BigInt::from(2).pow(k as u32) / to_bigint(p);
        let block_to_add = (m_bigint.bits() - 2 * NB as u64 + 1) / 2;
        let m = from_bigint::<U512>(&m_bigint);
        let mut x =
            server_key.extend_radix_with_trivial_zero_blocks_msb(x, NB + block_to_add as usize);
        let mut q = server_key.smart_scalar_mul_parallelized(&mut x, m);
        server_key.scalar_right_shift_assign_parallelized(&mut q, k as u64);
        server_key
            .sub_assign_parallelized(&mut x, &server_key.smart_scalar_mul_parallelized(&mut q, p));
        let len = x.blocks().len();
        server_key.trim_radix_blocks_msb_assign(&mut x, len - (NB + 1));

        return modulo_fast::<NB, _>(&x, p, server_key);
    }
    let c_blocks = (c.bits() as usize + 1) / 2;
    let x = server_key.extend_radix_with_trivial_zero_blocks_msb(x, (NB * 2) - x.blocks().len());

    // first pass NB*2 blocks
    let x_mod_p = (|x: &RadixCiphertext| {
        let mut a = server_key.scalar_right_shift_parallelized(x, n as u64);
        let mut b = server_key.smart_sub_parallelized(
            &mut x.clone(),
            &mut server_key.scalar_left_shift_parallelized(&a, n as u64),
        );

        let len = x.blocks().len();
        // a will be multiplied by c, so it must be at least NB + c_blocks long
        server_key.trim_radix_blocks_msb_assign(&mut a, len - (NB + c_blocks));
        // b must be at least NB long
        server_key.trim_radix_blocks_msb_assign(&mut b, len - NB);
        let ca = server_key.smart_scalar_mul_parallelized(&mut a, bigint_to_u128(&c));
        server_key.add_parallelized(&ca, &b)
    })(&x);

    // second pass % NB + c_blocks blocks
    let x_mod_p2 = (|x: &RadixCiphertext| {
        let mut a = server_key.scalar_right_shift_parallelized(x, n as u64);
        let mut b = server_key.smart_sub_parallelized(
            &mut x.clone(),
            &mut server_key.scalar_left_shift_parallelized(&a, n as u64),
        );

        let len = x.blocks().len();
        // a will be multiplied by c, so it must be at least NB + 1 long
        server_key.trim_radix_blocks_msb_assign(&mut a, len - (NB + 1));
        // b must be at least NB long
        server_key.trim_radix_blocks_msb_assign(&mut b, len - NB);
        let ca = server_key.smart_scalar_mul_parallelized(&mut a, bigint_to_u128(&c));
        server_key.add_parallelized(&ca, &b)
    })(&x_mod_p);

    modulo_fast::<NB, _>(&x_mod_p2, p, server_key)
}

/// Calculate a * b mod p
pub fn mul_mod_mersenne<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let mut a_expanded = server_key.extend_radix_with_trivial_zero_blocks_msb(a, NB);
    server_key.smart_mul_assign_parallelized(&mut a_expanded, &mut b.clone());
    mod_mersenne::<NB, _>(&a_expanded, p, server_key)
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use num_bigint::BigInt;
    use tfhe::{integer::keycache::IntegerKeyCache, shortint::prelude::PARAM_MESSAGE_2_CARRY_2};

    use crate::ops::{
        mersenne::{mersenne_mod_native, mul_mod_mersenne},
        native::mul_mod_native,
    };

    use super::{mersenne_coeff, mersenne_coeff_p};

    #[test]
    fn correct_mersenne_native_mod() {
        let p: u128 = 251;
        let x: u128 = 249;
        let y: u128 = 250;
        assert_eq!(mersenne_mod_native(x * y, p), mul_mod_native(x, y, p));

        let x = 123;
        let y = 250;
        assert_eq!(mersenne_mod_native(x * y, p), mul_mod_native(x, y, p));
    }

    #[test]
    fn correct_mersenne_mul_mod() {
        let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);
        const NUM_BLOCK: usize = 4;
        let p: u8 = 251;

        let mul_mod_naive = |x: u128, y: u128| -> u128 { (x * y) % p as u128 };

        let x: u128 = 250;
        let y: u128 = 249;
        let enc_x = client_key.encrypt_radix(x, NUM_BLOCK);
        let enc_y = client_key.encrypt_radix(y, NUM_BLOCK);
        let now = Instant::now();
        let xy_mod_p = mul_mod_mersenne::<NUM_BLOCK, _>(&enc_x, &enc_y, p, &server_key);
        println!(
            "mul mod mersenne done in {:.2}s",
            now.elapsed().as_secs_f64()
        );
        assert_eq!(
            client_key.decrypt_radix::<u128>(&xy_mod_p),
            mul_mod_naive(x, y)
        );
    }

    #[test]
    fn correct_mersenne_transfrom() {
        let p: u8 = 127;
        let coeff = mersenne_coeff_p(p);
        assert_eq!(coeff, (7, BigInt::from(1)));
    }
}
