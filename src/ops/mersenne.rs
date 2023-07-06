use num_bigint::BigInt;
use tfhe::integer::{block_decomposition::DecomposableInto, RadixCiphertext, ServerKey};

use crate::helper::bigint_to_u256;

/// Calculate n, m, p from coeff
#[inline(always)]
pub fn mersenne_coeff(coeff: &[u32]) -> (u32, Option<u32>, BigInt, BigInt, BigInt) {
    assert!(coeff.len() > 1);
    let len = coeff.len();
    let n = coeff[0];
    let m = match len > 2 {
        true => Some(coeff[1]),
        false => None,
    };
    let p = coeff[1..len - 1]
        .iter()
        .fold(BigInt::from(2).pow(coeff[0]), |acc, b| {
            acc - BigInt::from(2).pow(*b)
        })
        - coeff[len - 1];
    let q = BigInt::from(2).pow(n);
    let c = &p - BigInt::from(2).pow(n)
        + match m {
            Some(m) => BigInt::from(2).pow(m),
            None => BigInt::from(0),
        };

    (n, m, p, q, c)
}

/// Calculate x mod q mod p where q = 2^n_0
/// `coeff` in the form of p = 2^n_0 - 2^n_1 - ... - 2^n_{k-1} - n_k
pub fn mersenne_mod_native(coeff: &[u32], x: &BigInt) -> BigInt {
    assert!(coeff.len() > 1);
    let len = coeff.len();
    let n = coeff[0];
    let m = match len > 2 {
        true => Some(coeff[1]),
        false => None,
    };
    let p = coeff[1..len - 1]
        .iter()
        .fold(BigInt::from(2).pow(coeff[0]), |acc, b| {
            println!("acc: {acc}, b: {b}");
            acc - BigInt::from(2).pow(*b)
        })
        - coeff[len - 1];
    let q = BigInt::from(2).pow(n);
    let c = q + match m {
        Some(m) => BigInt::from(2).pow(m),
        None => BigInt::from(0),
    } - &p;

    // x = a*2^n + b
    let a = x >> n;
    let b = x - (&a << n);
    assert_eq!(*x, &a * &BigInt::from(2).pow(n) + &b);

    // x % p = (a*2^n + b) % p = a * 2^m + c * a + b % p
    let x_mod_p = match m {
        Some(m) => &a * &BigInt::from(2).pow(m),
        None => BigInt::from(0),
    } + (&a * &c)
        + (&b);
    assert_eq!(x_mod_p, x % &p);

    x_mod_p
}

pub fn mersenne_mod<const NB: usize>(
    x: &RadixCiphertext,
    p_coeff: &[u32],
    server_key: &ServerKey,
) -> RadixCiphertext {
    assert!(NB * 2 >= p_coeff[0] as usize);
    let (n, m, _, _, c) = mersenne_coeff(p_coeff);

    let a = server_key.scalar_right_shift_parallelized(&x, n as u64);
    let b =
        server_key.sub_parallelized(&x, &server_key.scalar_left_shift_parallelized(&a, n as u64));
    let a2m = m.map(|m| server_key.scalar_left_shift_parallelized(&a, m as u64));
    let cab = server_key.add_parallelized(
        &server_key.mul_parallelized(&a, &server_key.create_trivial_radix(bigint_to_u256(&c), NB)),
        &b,
    );

    match a2m {
        Some(a2m) => server_key.add_parallelized(&a2m, &cab),
        None => cab,
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use tfhe::{
        integer::{keycache::IntegerKeyCache, U256},
        shortint::prelude::PARAM_MESSAGE_2_CARRY_2,
    };

    use crate::{
        helper::{bigint_to_u256, u256_to_bigint},
        ops::mersenne::{mersenne_mod, mersenne_mod_native},
    };

    use super::mersenne_coeff;

    #[test]
    fn correct_mersenne_native_mod() {
        let coeff: &[u32] = &[7, 1];
        let (_, _, p, q, _) = mersenne_coeff(coeff);
        let x = BigInt::from(253) % &q;

        assert_eq!(mersenne_mod_native(coeff, &x), &x % &p);
    }

    #[test]
    fn correct_mersenne_mod() {
        let coeff: &[u32] = &[7, 1];
        let (_, _, p, q, _) = mersenne_coeff(coeff);
        let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);
        const NUM_BLOCK: usize = 4;

        let x = BigInt::from(253) % &q;
        let enc_x = client_key.encrypt_radix(bigint_to_u256(&x), NUM_BLOCK);
        let enc_x_mod_p = mersenne_mod::<NUM_BLOCK>(&enc_x, coeff, &server_key);
        let dec_x_mod_p = u256_to_bigint(client_key.decrypt_radix::<U256>(&enc_x_mod_p));
        assert_eq!(dec_x_mod_p, &x % &p);

        let y = BigInt::from(157) % &q;
        let enc_y = client_key.encrypt_radix(bigint_to_u256(&y), NUM_BLOCK);
        let enc_y_mod_p = mersenne_mod::<NUM_BLOCK>(&enc_y, coeff, &server_key);
        let dec_y_mod_p = u256_to_bigint(client_key.decrypt_radix::<U256>(&enc_y_mod_p));
        assert_eq!(dec_y_mod_p, &y % &p);

        let z = BigInt::from(0) % &q;
        let enc_z = client_key.encrypt_radix(bigint_to_u256(&z), NUM_BLOCK);
        let enc_z_mod_p = mersenne_mod::<NUM_BLOCK>(&enc_z, coeff, &server_key);
        let dec_z_mod_p = u256_to_bigint(client_key.decrypt_radix::<U256>(&enc_z_mod_p));
        assert_eq!(dec_z_mod_p, &z % &p);
    }
}
