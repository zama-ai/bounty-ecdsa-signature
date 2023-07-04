use std::time::Instant;

use tfhe::{
    core_crypto::prelude::{Numeric, UnsignedInteger},
    integer::ClientKey,
    integer::{
        block_decomposition::{DecomposableInto, RecomposableFrom},
        IntegerCiphertext, RadixCiphertext, ServerKey,
    },
};

use crate::helper::format;

pub mod group;

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
    for i in 0..<P as Numeric>::BITS {
        let now = Instant::now();
        // q, r = r0 / r1
        let (mut q, mut r) =
            server_key.smart_div_rem_parallelized(&mut r0.clone(), &mut r1.clone());
        rayon::join(
            || server_key.full_propagate_parallelized(&mut q),
            || server_key.full_propagate_parallelized(&mut r),
        );
        let tmp = t1.clone();
        // t1 = t0 - q * t1
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
        // inv = inv + done_now * t1
        let mut update = server_key.smart_mul_parallelized(&mut done_now, &mut t0);
        server_key.smart_add_assign_parallelized(&mut inv, &mut update);
        // update values
        r0 = r1;
        r1 = r;

        println!(
            "Inverse mod bit {i} took {:.2}s",
            now.elapsed().as_secs_f64()
        );
    }

    // final result mod p
    // inverse can be **negative**. so we need to add p to make it positive
    server_key.extend_radix_with_trivial_zero_blocks_msb_assign(&mut inv, 1);
    server_key.smart_scalar_add_assign_parallelized(&mut inv, p);
    let mut is_gt = server_key.smart_scalar_ge_parallelized(&mut inv, p);
    server_key.trim_radix_blocks_msb_assign(&mut is_gt, NB - 1);
    let mut to_sub =
        server_key.smart_mul_parallelized(&mut server_key.create_trivial_radix(p, NB), &mut is_gt);
    server_key.smart_sub_assign_parallelized(&mut inv, &mut to_sub);
    server_key.full_propagate_parallelized(&mut inv);
    server_key.trim_radix_blocks_msb_assign(&mut inv, 1);
    inv
}

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
    server_key.full_propagate_parallelized(&mut a_expanded);
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
    server_key.full_propagate_parallelized(&mut a_expanded);
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

        println!("inverse bit {i} took {:.2}s", now.elapsed().as_secs_f32());
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
    server_key.full_propagate_parallelized(&mut a_expanded);
    let (_q, mut r) = server_key.smart_div_rem_parallelized(
        &mut a_expanded,
        &mut server_key.create_trivial_radix(p, NB * 2),
    );
    server_key.full_propagate_parallelized(&mut r);
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
        ops::{add_mod, double_mod, mul_mod, mul_mod_constant, sub_mod},
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
    fn correct_add_mod() {
        let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);
        const NUM_BLOCK: usize = 4;
        let p: u8 = 251;

        let a: u128 = 248;
        let b: u128 = 249;
        let c: u128 = (a + b) % p as u128;
        let enc_c = add_mod::<NUM_BLOCK, _>(
            &client_key.encrypt_radix(a, NUM_BLOCK),
            &client_key.encrypt_radix(b, NUM_BLOCK),
            p,
            &server_key,
        );
        assert_eq!(c as u8, client_key.decrypt_radix::<u8>(&enc_c));

        let d = (c + c) % p as u128;
        let enc_d = add_mod::<NUM_BLOCK, _>(&enc_c, &enc_c, p, &server_key);
        assert_eq!(d as u8, client_key.decrypt_radix::<u8>(&enc_d));

        let e = (c + a) % p as u128;
        let enc_e = add_mod::<NUM_BLOCK, _>(
            &enc_c,
            &client_key.encrypt_radix(a, NUM_BLOCK),
            p,
            &server_key,
        );
        assert_eq!(e as u8, client_key.decrypt_radix::<u8>(&enc_e));

        let f = (e + b) % p as u128;
        let enc_f = add_mod::<NUM_BLOCK, _>(
            &enc_e,
            &client_key.encrypt_radix(b, NUM_BLOCK),
            p,
            &server_key,
        );
        assert_eq!(f as u8, client_key.decrypt_radix::<u8>(&enc_f));
    }

    #[test]
    fn correct_sub_mod() {
        let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);
        const NUM_BLOCK: usize = 4;
        let p: u8 = 251;

        let sub_mod_naive = |a: u128, b: u128| (a + p as u128 - b) % p as u128;

        let a: u128 = 248;
        let b: u128 = 249;
        let c: u128 = sub_mod_naive(a, b);
        let enc_c = sub_mod::<NUM_BLOCK, _>(
            &client_key.encrypt_radix(a, NUM_BLOCK),
            &client_key.encrypt_radix(b, NUM_BLOCK),
            p,
            &server_key,
        );
        assert_eq!(c as u8, client_key.decrypt_radix::<u8>(&enc_c));

        let d = sub_mod_naive(c, b);
        let enc_d = sub_mod::<NUM_BLOCK, _>(
            &enc_c,
            &client_key.encrypt_radix(b, NUM_BLOCK),
            p,
            &server_key,
        );
        assert_eq!(d as u8, client_key.decrypt_radix::<u8>(&enc_d));

        let e = sub_mod_naive(c, d);
        let enc_e = sub_mod::<NUM_BLOCK, _>(&enc_c, &enc_d, p, &server_key);
        assert_eq!(e as u8, client_key.decrypt_radix::<u8>(&enc_e));
    }

    #[test]
    fn correct_mul_mod() {
        let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);
        const NUM_BLOCK: usize = 4;
        let p: u8 = 251;

        let mul_mod_naive = |a: u128, b: u128| (a * b) % p as u128;

        let a: u128 = 248;
        let b: u128 = 249;
        let c: u128 = mul_mod_naive(a, b);
        let enc_c = mul_mod::<NUM_BLOCK, _>(
            &client_key.encrypt_radix(a, NUM_BLOCK),
            &client_key.encrypt_radix(b, NUM_BLOCK),
            p,
            &server_key,
        );
        assert_eq!(c as u8, client_key.decrypt_radix::<u8>(&enc_c));

        let d = mul_mod_naive(c, b);
        let enc_d = mul_mod::<NUM_BLOCK, _>(
            &enc_c,
            &client_key.encrypt_radix(b, NUM_BLOCK),
            p,
            &server_key,
        );
        assert_eq!(d as u8, client_key.decrypt_radix::<u8>(&enc_d));

        let e = mul_mod_naive(c, d);
        let enc_e = mul_mod::<NUM_BLOCK, _>(&enc_c, &enc_d, p, &server_key);
        assert_eq!(e as u8, client_key.decrypt_radix::<u8>(&enc_e));

        let f = mul_mod_naive(e, e);
        let enc_f = mul_mod::<NUM_BLOCK, _>(&enc_e, &enc_e, p, &server_key);
        assert_eq!(f as u8, client_key.decrypt_radix::<u8>(&enc_f));

        let g = mul_mod_naive(f, f);
        let enc_g = mul_mod::<NUM_BLOCK, _>(&enc_f, &enc_f, p, &server_key);
        assert_eq!(g as u8, client_key.decrypt_radix::<u8>(&enc_g));
    }
}
