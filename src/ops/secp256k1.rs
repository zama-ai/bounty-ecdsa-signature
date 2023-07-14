use std::time::Instant;

use num_bigint::BigInt;
use rand::Rng;
use tfhe::integer::{
    block_decomposition::{DecomposableInto, RecomposableFrom},
    IntegerCiphertext, RadixCiphertext, ServerKey,
};

use crate::ops::{add_mod, double_mod, mul_mod, square_mod, sub_mod};

/// Calculate x mod p^2 mod p
/// where p is the prime modulus of secp256k1
pub fn modulo_native(x: &BigInt, p: &BigInt) -> BigInt {
    let set = |x: &BigInt, n: usize| {
        let x = x % &BigInt::from(2).pow(256 - n as u32);
        x << n
    };
    let c0 = x % &BigInt::from(2).pow(256);
    let c1: BigInt = x >> 256;
    let w1 = set(&c1, 32);
    let w2 = set(&c1, 9);
    let w3 = set(&c1, 8);
    let w4 = set(&c1, 7);
    let w5 = set(&c1, 6);
    let w6 = set(&c1, 4);
    let w7 = c1.clone();

    let k1 = (&c1 >> 252) + (&c1 >> 250);
    let k2 = &k1 + (&c1 >> 249);
    let k3 = &k2 + (&c1 >> 248);
    let k4 = &k3 + (&c1 >> 247);
    let s1 = &k4 + (&c1 >> 224);
    let k11 = ((&s1 << 2) + (&s1 << 1)) + &s1;
    let k12 = &k11 << 7;
    let k13 = &(&s1 << 4) + &s1;
    let k14 = &(&s1 << 6) + &k13;
    let k = &(&s1 << 32) + &k12 + &k14;
    let s = &c0 + &k + &w1 + &w2 + &w3 + &w4 + &w5 + &w6 + &w7;

    s
    //match &s > p {
    //true => s - p,
    //false => s,
    //}
}

/// Calculate x mod p^2 mod p
/// where p is the prime modulus of secp256k1
pub fn modulo<const NB: usize>(x: &RadixCiphertext, server_key: &ServerKey) -> RadixCiphertext {
    #[cfg(feature = "low_level_timing")]
    let ops_start = Instant::now();
    #[cfg(feature = "low_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);

    #[cfg(feature = "low_level_timing")]
    println!("kem modulo start -- ref {}", task_ref);

    let mut x = x.clone();
    let len = x.blocks().len();
    if len < 2 * NB {
        server_key.extend_radix_with_trivial_zero_blocks_msb_assign(&mut x, 2 * NB - len);
    }

    let c0 = RadixCiphertext::from_blocks(x.blocks()[..NB].to_vec());
    let c1 = RadixCiphertext::from_blocks(x.blocks()[NB..].to_vec());

    let w1 = server_key.scalar_left_shift_parallelized(&c1, 32);
    let w2 = server_key.scalar_left_shift_parallelized(&c1, 9);
    let w3 = server_key.scalar_left_shift_parallelized(&c1, 8);
    let w4 = server_key.scalar_left_shift_parallelized(&c1, 7);
    let w5 = server_key.scalar_left_shift_parallelized(&c1, 6);
    let w6 = server_key.scalar_left_shift_parallelized(&c1, 4);
    let w7 = c1.clone();

    let k1 = {
        let lhs = server_key.scalar_right_shift_parallelized(&c1, 252);
        let rhs = server_key.scalar_right_shift_parallelized(&c1, 250);
        server_key.add_parallelized(&lhs, &rhs)
    };
    let k2 = {
        let rhs = server_key.scalar_right_shift_parallelized(&c1, 249);
        server_key.add_parallelized(&k1, &rhs)
    };
    let k3 = {
        let rhs = server_key.scalar_right_shift_parallelized(&c1, 248);
        server_key.add_parallelized(&k2, &rhs)
    };
    let k4 = {
        let rhs = server_key.scalar_right_shift_parallelized(&c1, 247);
        server_key.add_parallelized(&k3, &rhs)
    };
    let s1 = {
        let rhs = server_key.scalar_right_shift_parallelized(&c1, 224);
        server_key.add_parallelized(&k4, &rhs)
    };
    let k11 = {
        let rhs = server_key.scalar_left_shift_parallelized(&s1, 2);
        let lhs = server_key.scalar_left_shift_parallelized(&s1, 1);
        server_key.add_parallelized(&server_key.add_parallelized(&rhs, &lhs), &s1)
    };
    let k12 = server_key.scalar_left_shift_parallelized(&k11, 7);
    let k13 = server_key.add_parallelized(&server_key.scalar_left_shift_parallelized(&s1, 4), &s1);
    let k14 = server_key.add_parallelized(&server_key.scalar_left_shift_parallelized(&s1, 6), &k13);
    let k = server_key.add_parallelized(
        &server_key.scalar_left_shift_parallelized(&s1, 32),
        &server_key.add_parallelized(&k12, &k14),
    );
    let s = {
        let t0 = server_key.add_parallelized(&c0, &k);
        let t1 = server_key.add_parallelized(&t0, &w1);
        let t2 = server_key.add_parallelized(&t1, &w2);
        let t3 = server_key.add_parallelized(&t2, &w3);
        let t4 = server_key.add_parallelized(&t3, &w4);
        let t5 = server_key.add_parallelized(&t4, &w5);
        let t6 = server_key.add_parallelized(&t5, &w6);
        server_key.add_parallelized(&t6, &w7)
    };

    #[cfg(feature = "low_level_timing")]
    println!(
        "kem modulo done in {:.2}s -- ref {}",
        ops_start.elapsed().as_secs_f64(),
        task_ref
    );
    s
}

#[cfg(test)]
mod tests {
    use std::{str::FromStr, time::Instant};

    use num_bigint::BigInt;
    use tfhe::{
        integer::{keycache::IntegerKeyCache, U256},
        shortint::prelude::PARAM_MESSAGE_2_CARRY_2,
    };

    use crate::{
        helper::{format, from_bigint, to_bigint, u256_from_decimal_string},
        ops::{group_jacobian::group_projective_into_affine, secp256k1::modulo},
    };

    use super::modulo_native;

    #[test]
    fn correct_secp_reduction_native() {
        let p = BigInt::from_str(
            "115792089237316195423570985008687907853269984665640564039457584007908834671663",
        )
        .unwrap();
        let x1 = BigInt::from_str(
            "55066263022277343669578718895168534326250603453777594175500187360389116729240",
        )
        .unwrap();
        let y1 = BigInt::from_str(
            "32670510020758816978083085130507043184471273380659243275938904335757337482424",
        )
        .unwrap();

        let mul_mod_naive = |x: &BigInt, y: &BigInt| {
            let mut res = x * y;
            res %= &p;
            res
        };
        let mul_mod = |x: &BigInt, y: &BigInt| {
            let res = x * y;
            modulo_native(&res, &p)
        };

        let result = mul_mod(&x1, &y1);
        let result_naive = mul_mod_naive(&x1, &y1);
        assert_eq!(result, result_naive);
    }

    #[test]
    fn correct_secp_reduction() {
        let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);

        const NUM_BLOCK: usize = 128;
        type Integer = U256;
        let p: Integer = u256_from_decimal_string(
            "115792089237316195423570985008687907853269984665640564039457584007908834671663",
        );
        let x1: Integer = u256_from_decimal_string(
            "55066263022277343669578718895168534326250603453777594175500187360389116729240",
        );
        let y1: Integer = u256_from_decimal_string(
            "32670510020758816978083085130507043184471273380659243275938904335757337482424",
        );

        let ct_x1 = client_key.encrypt_radix(x1, NUM_BLOCK);
        let mut ct_y1 = client_key.encrypt_radix(y1, NUM_BLOCK);

        let mul_mod_int = |x: Integer, y: Integer| -> Integer {
            let x = to_bigint(x);
            let y = to_bigint(y);
            let p = to_bigint(p);
            from_bigint(&((&x * &y) % &p))
        };

        let now = Instant::now();
        let mut a_expanded =
            server_key.extend_radix_with_trivial_zero_blocks_msb(&ct_x1, NUM_BLOCK);
        server_key.smart_mul_assign_parallelized(&mut a_expanded, &mut ct_y1);
        let res = modulo::<NUM_BLOCK>(&a_expanded, &server_key);
        println!(
            "mul and secp reduction done in {:.2}s",
            now.elapsed().as_secs_f64()
        );
        println!("got: {}", format(client_key.decrypt_radix::<U256>(&res)));
        println!("should be: {}", format(mul_mod_int(x1, y1)));
    }
}
