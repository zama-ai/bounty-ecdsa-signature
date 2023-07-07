use num_bigint::BigInt;
use tfhe::integer::{
    block_decomposition::{DecomposableInto, RecomposableFrom},
    IntegerCiphertext, RadixCiphertext, ServerKey, U256,
};

use crate::helper::{
    bigint_ilog2_ceil, bigint_to_u256, format, read_client_key, to_bigint, u256_to_bigint,
};

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
pub fn mersenne_coeff_p<P: DecomposableInto<u8> + Copy + Sync>(p: P) -> (u32, BigInt) {
    let pb = to_bigint(p);
    let n = bigint_ilog2_ceil(&pb);
    let c = (BigInt::from(1) << n) - &pb;

    (n, c)
}

/// Calculate x mod p^2 mod p
/// `coeff` in the form of p = 2^n_0 - 2^n_1 - ... - 2^n_{k-1} - n_k
pub fn mersenne_mod_native(coeff: &[u32], x: &BigInt) -> BigInt {
    let (n, p, _, c) = mersenne_coeff(coeff);

    // x = a*2^n + b
    let a = x >> n;
    let b = x - (&a << n);
    assert_eq!(*x, &a * &BigInt::from(2).pow(n) + &b);

    // x % p = (a*2^n + b) % p = c * a + b % p
    let x_mod_p = &a * &c + &b;
    assert_eq!(x_mod_p, x % &p);

    x_mod_p
}

/// Calculate x mod p^2 mod p
/// `coeff` in the form of p = 2^n_0 - 2^n_1 - ... - 2^n_{k-1} - n_k
pub fn mod_mersenne_coeff<const NB: usize>(
    x: &RadixCiphertext,
    p_coeff: &[u32],
    server_key: &ServerKey,
) -> RadixCiphertext {
    let (n, p, _, c) = mersenne_coeff(p_coeff);

    let mut a = server_key.scalar_right_shift_parallelized(&x, n as u64);
    let mut b =
        server_key.sub_parallelized(&x, &server_key.scalar_left_shift_parallelized(&a, n as u64));
    let len = a.blocks().len();
    if len > NB {
        server_key.trim_radix_blocks_msb_assign(&mut a, len - NB);
        server_key.trim_radix_blocks_msb_assign(&mut b, len - NB);
    }
    // c * a + b
    let mut cab = server_key.add_parallelized(
        &b,
        &server_key.mul_parallelized(&a, &server_key.create_trivial_radix(bigint_to_u256(&c), NB)),
    );
    let len = cab.blocks().len();
    server_key.trim_radix_blocks_msb_assign(&mut cab, len - NB);
    // if cab >= p, cab -= p
    let mut is_gt = server_key.scalar_gt_parallelized(&cab, bigint_to_u256(&p));
    server_key.trim_radix_blocks_msb_assign(&mut is_gt, NB - 1);
    server_key.smart_sub_assign_parallelized(
        &mut cab,
        &mut server_key.smart_mul_parallelized(
            &mut server_key.create_trivial_radix(bigint_to_u256(&p), NB),
            &mut is_gt,
        ),
    );
    server_key.full_propagate_parallelized(&mut cab);
    cab
}

/// Calculate x mod p^2 mod p
pub fn mod_mersenne<
    const NB: usize,
    P: DecomposableInto<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    x: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let (n, c) = mersenne_coeff_p(p);
    let process = |x: &RadixCiphertext| {
        let mut a = server_key.scalar_right_shift_parallelized(x, n as u64);
        let len = x.blocks().len();
        let mut b = server_key.smart_sub_parallelized(
            &mut x.clone(),
            &mut server_key.scalar_left_shift_parallelized(&a, n as u64),
        );
        server_key.trim_radix_blocks_msb_assign(&mut a, len - NB);
        server_key.trim_radix_blocks_msb_assign(&mut b, len - NB);
        let ca = server_key.mul_parallelized(
            &server_key.create_trivial_radix(bigint_to_u256(&c), NB * 3 / 2),
            &a,
        );
        let x_mod_p = server_key.add_parallelized(&ca, &b);
        x_mod_p
    };
    let x_mod_p = process(x);
    let x_mod_p2 = process(&x_mod_p);
    let mut x_mod_p3 = process(&x_mod_p2);
    let len = x_mod_p3.blocks().len();
    server_key.trim_radix_blocks_msb_assign(&mut x_mod_p3, len - NB);
    x_mod_p3
}

/// Calculate a * b mod p
/// `coeff` in the form of p = 2^n_0 - 2^n_1 - ... - 2^n_{k-1} - n_k
pub fn mul_mod_mersenne_coeff<const NB: usize>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p_coeff: &[u32],
    server_key: &ServerKey,
) -> RadixCiphertext {
    let mut a_expanded = server_key.extend_radix_with_trivial_zero_blocks_msb(a, NB);
    server_key.smart_mul_assign_parallelized(&mut a_expanded, &mut b.clone());
    server_key.full_propagate_parallelized(&mut a_expanded);
    mod_mersenne_coeff::<NB>(&a_expanded, p_coeff, server_key)
}

/// Calculate a * b mod p
pub fn mul_mod_mersenne<
    const NB: usize,
    P: DecomposableInto<u64> + RecomposableFrom<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let mut a_expanded = server_key.extend_radix_with_trivial_zero_blocks_msb(a, NB);
    server_key.smart_mul_assign_parallelized(&mut a_expanded, &mut b.clone());
    server_key.full_propagate_parallelized(&mut a_expanded);
    mod_mersenne::<NB, _>(&a_expanded, p, server_key)
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use tfhe::{
        integer::{keycache::IntegerKeyCache, U256},
        shortint::prelude::PARAM_MESSAGE_2_CARRY_2,
    };

    use crate::{
        helper::{bigint_to_u256, set_client_key, u256_to_bigint},
        ops::mersenne::{mersenne_mod_native, mod_mersenne_coeff, mul_mod_mersenne_coeff},
    };

    use super::{mersenne_coeff, mersenne_coeff_p};

    #[test]
    fn correct_mersenne_native_mod() {
        let coeff: &[u32] = &[7, 1];
        let (_, p, _q, _) = mersenne_coeff(coeff);
        let p2 = p.pow(2);
        let x = BigInt::from(16000) % &p2;

        assert_eq!(mersenne_mod_native(coeff, &x), &x % &p);
    }

    #[test]
    fn correct_mersenne_mod() {
        let coeff: &[u32] = &[7, 1];
        let (_, p, _q, _) = mersenne_coeff(coeff);
        let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);
        const NUM_BLOCK: usize = 7;
        let p2 = p.pow(2);

        let x = BigInt::from(16000) % &p2;
        let enc_x = client_key.encrypt_radix(bigint_to_u256(&x), NUM_BLOCK);
        let enc_x_mod_p = mod_mersenne_coeff::<NUM_BLOCK>(&enc_x, coeff, &server_key);
        let dec_x_mod_p = u256_to_bigint(client_key.decrypt_radix::<U256>(&enc_x_mod_p));
        assert_eq!(dec_x_mod_p, &x % &p);

        let y = BigInt::from(1500) % &p2;
        let enc_y = client_key.encrypt_radix(bigint_to_u256(&y), NUM_BLOCK);
        let enc_y_mod_p = mod_mersenne_coeff::<NUM_BLOCK>(&enc_y, coeff, &server_key);
        let dec_y_mod_p = u256_to_bigint(client_key.decrypt_radix::<U256>(&enc_y_mod_p));
        assert_eq!(dec_y_mod_p, &y % &p);

        let z = BigInt::from(0) % &p2;
        let enc_z = client_key.encrypt_radix(bigint_to_u256(&z), NUM_BLOCK);
        let enc_z_mod_p = mod_mersenne_coeff::<NUM_BLOCK>(&enc_z, coeff, &server_key);
        let dec_z_mod_p = u256_to_bigint(client_key.decrypt_radix::<U256>(&enc_z_mod_p));
        assert_eq!(dec_z_mod_p, &z % &p);
    }

    #[test]
    fn correct_mersenne_mul_mod() {
        let coeff: &[u32] = &[7, 1];
        let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);
        set_client_key(&client_key);
        const NUM_BLOCK: usize = 4;
        let p: u8 = 127;

        let mul_mod_naive = |x: u128, y: u128| -> u128 { (x * y) % p as u128 };

        let x = 125;
        let y = 126;
        let enc_x = client_key.encrypt_radix(x, NUM_BLOCK);
        let enc_y = client_key.encrypt_radix(y, NUM_BLOCK);
        let xy_mod_p = mul_mod_mersenne_coeff::<NUM_BLOCK>(&enc_x, &enc_y, coeff, &server_key);
        assert_eq!(
            client_key.decrypt_radix::<u128>(&xy_mod_p),
            mul_mod_naive(x, y)
        );
    }

    #[test]
    fn correct_mersenne_transfrom() {
        let p = 127;
        let coeff = mersenne_coeff_p(p);
        assert_eq!(coeff, (7, BigInt::from(1)));
    }
}
