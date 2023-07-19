use std::time::Instant;

use rand::Rng;
use tfhe::{
    core_crypto::prelude::{Numeric, UnsignedInteger},
    integer::{
        block_decomposition::{DecomposableInto, RecomposableFrom},
        IntegerCiphertext, RadixCiphertext, ServerKey,
    },
};

use crate::{
    helper::{format, read_client_key},
    numeral::Numeral,
    ops::mersenne::{mod_mersenne, mod_mersenne_fast},
    stats::{ProtocolLowOps, ProtocolStats},
};

use self::mersenne::mul_mod_mersenne;

pub mod group_homogenous;
pub mod group_jacobian;
pub mod mersenne;
pub mod native;
pub mod secp256k1;

/// selector ? a : 0
/// selector is a bit (0 or 1)
pub fn selector_zero(
    a: &RadixCiphertext,
    selector: &RadixCiphertext,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let mut res = a.clone();
    let len = selector.blocks().len();
    let mut selector = server_key.trim_radix_blocks_msb(&selector, len - 1);
    //let mut not_is_gt =
    //server_key.smart_sub_parallelized(&mut server_key.create_trivial_radix(0, len), &mut is_gt);
    server_key.smart_mul_assign_parallelized(&mut res, &mut selector);
    res
}

/// selector ? a : 0
/// selector is a bit (0 or 1)
pub fn selector_zero_constant<const NB: usize, P: Numeral>(
    a: P,
    selector: &RadixCiphertext,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let len = selector.blocks().len();
    #[cfg(not(feature = "inw_selector"))]
    let res = {
        let mut selector =
            server_key.extend_radix_with_trivial_zero_blocks_msb(&selector, NB - len);
        server_key.smart_scalar_mul_assign_parallelized(&mut selector, a);
        selector
    };
    #[cfg(feature = "inw_selector")]
    let res = {
        let mut not_selector = server_key.smart_sub_parallelized(
            &mut server_key.create_trivial_radix(0, NB),
            &mut server_key.extend_radix_with_trivial_zero_blocks_msb(&selector, NB - len),
        );
        server_key.smart_scalar_bitand_assign_parallelized(&mut not_selector, a);
        not_selector
    };
    res
}

/// selector ? a : b
/// selector is a bit (0 or 1)
pub fn selector(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    selector: &RadixCiphertext,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let len = selector.blocks().len();
    let mut selector = server_key.trim_radix_blocks_msb(&selector, len - 1);
    let mut not_selector = server_key
        .smart_sub_parallelized(&mut server_key.create_trivial_radix(1, len), &mut selector);
    let (mut r0, mut r1) = rayon::join(
        || server_key.smart_mul_parallelized(&mut a.clone(), &mut selector),
        || server_key.smart_mul_parallelized(&mut b.clone(), &mut not_selector),
    );
    server_key.smart_add_parallelized(&mut r0, &mut r1)
}

/// a_0 + a_1 + ... + a_n mod p
pub fn multi_add_mod<const NB: usize, P: Numeral>(
    a: &[RadixCiphertext],
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    // assume large p and a,b < p

    #[cfg(feature = "low_level_timing")]
    let start_ops = Instant::now();
    #[cfg(feature = "low_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);

    // add all elements in a together
    let le = a.len();
    // find bit length of le
    let log_le = le.ceil_ilog2() as usize;
    // find block length of le
    let extend = log_le / 2 + 1;

    // create sum object
    let mut sum = server_key.extend_radix_with_trivial_zero_blocks_msb(&a[0], extend);

    // just sum all elements
    for ai in a {
        //let mut tmp = server_key.extend_radix_with_trivial_zero_blocks_msb(&a[i].clone(), extend);
        server_key.smart_add_assign_parallelized(&mut sum, &mut ai.clone());
    }
    // only check with p*2^i from high to low
    for i in (0..le).rev() {
        // to_check = p * 2^i = p << i
        let mut to_check = server_key.create_trivial_radix(p, NB + extend);
        server_key.scalar_left_shift_assign_parallelized(&mut to_check, i as u64);
        let mut is_gt = server_key.smart_gt_parallelized(&mut sum, &mut to_check);
        server_key.trim_radix_blocks_msb_assign(&mut is_gt, NB + extend - 1);
        let mut to_sub = server_key.smart_mul_parallelized(&mut to_check, &mut is_gt);
        server_key.smart_sub_assign_parallelized(&mut sum, &mut to_sub);
    }

    // server_key.full_propagate_parallelized(&mut sum);
    #[cfg(feature = "low_level_timing")]
    println!(
        "multi add mod done in {:.2}s -- ref {}",
        start_ops.elapsed().as_secs_f64(),
        task_ref
    );

    sum
}

/// turn x mod a to x mod b
/// only if a > b and a < 2b
pub fn modulo_fast<const NB: usize, P: Numeral>(
    x: &RadixCiphertext,
    b: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let len = x.blocks().len();
    let mut x = x.clone();
    let is_gt = server_key.smart_scalar_ge_parallelized(&mut x, b);
    let to_sub = selector_zero_constant::<NB, _>(b, &is_gt, server_key);
    server_key.sub_assign_parallelized(&mut x, &to_sub);
    server_key.trim_radix_blocks_msb_assign(&mut x, len - NB);
    x
}

/// turn x mod a to x mod b
/// for all cases, require 1 division
pub fn modulo_div_rem<
    const NB: usize,
    P: DecomposableInto<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    x: &RadixCiphertext,
    b: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let (_q, r) = server_key
        .smart_div_rem_parallelized(&mut x.clone(), &mut server_key.create_trivial_radix(b, NB));
    // server_key.full_propagate_parallelized(&mut r);
    r
}

pub fn inverse_mod_binary_gcd<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    // implement binary gcd algorithm
    // assume a < p. (no check)
    // This is a draft algo. TODO: remove hard coded test value.
    let mut a = a.clone();
    let mut u: tfhe::integer::ciphertext::BaseRadixCiphertext<tfhe::shortint::Ciphertext> =
        server_key.create_trivial_radix(1, NB);
    let mut b = server_key.create_trivial_radix(p, NB);
    let mut v = server_key.create_trivial_radix(0, NB);
    let mut was_done = server_key.create_trivial_radix(0, 1);
    let loop_end = <P as Numeric>::BITS * 2;
    let mut result = server_key.create_trivial_radix(0, NB);
    read_client_key(|server_key| {
        println!(
            "a = {}\nb = {}\nu = {}\nv = {}\n-----",
            format(server_key.decrypt_radix::<P>(&a)),
            format(server_key.decrypt_radix::<P>(&b)),
            format(server_key.decrypt_radix::<P>(&u)),
            format(server_key.decrypt_radix::<P>(&v)),
        )
    });

    for _ in 0..loop_end {
        // if a is even -> a = a/2, u = u/2
        let mut a_is_odd = server_key.scalar_bitand_parallelized(&a, 1);
        server_key.trim_radix_blocks_msb_assign(&mut a_is_odd, NB - 1);
        let a_is_even = server_key
            .smart_sub_parallelized(&mut server_key.create_trivial_radix(1, 1), &mut a_is_odd);

        let a_c1 = server_key.scalar_right_shift_parallelized(&a, 1);
        let b_c1 = b.clone();

        let u_div_2 = mul_mod::<NB, _>(
            &u,
            &server_key.create_trivial_radix(126u64, NB), // TODO: fix this
            p,
            server_key,
        ); // test for 251
        let u_c1 = mersenne::mod_mersenne::<NB, _>(&u_div_2, p, server_key);
        let v_c1 = v.clone();

        // if a < b then (a, u, b, v) ← (b, v, a, u)
        let mut a_lt_b = server_key.smart_lt_parallelized(&mut a.clone(), &mut b.clone());
        server_key.trim_radix_blocks_msb_assign(&mut a_lt_b, NB - 1);
        let mut a_ge_b = server_key
            .smart_sub_parallelized(&mut server_key.create_trivial_radix(1, 1), &mut a_lt_b);
        let (mut sa, mut su, mut sb, mut sv) = (b.clone(), v.clone(), a.clone(), u.clone());
        let sa_c2 = server_key.smart_add_parallelized(
            &mut server_key
                .smart_mul_parallelized(&mut sa, &mut a_lt_b)
                .clone(),
            &mut server_key
                .smart_mul_parallelized(&mut a, &mut a_ge_b)
                .clone(),
        );
        let sb_c2 = server_key.smart_add_parallelized(
            &mut server_key
                .smart_mul_parallelized(&mut sb, &mut a_lt_b)
                .clone(),
            &mut server_key
                .smart_mul_parallelized(&mut b, &mut a_ge_b)
                .clone(),
        );
        let su_c2 = server_key.smart_add_parallelized(
            &mut server_key
                .smart_mul_parallelized(&mut su, &mut a_lt_b)
                .clone(),
            &mut server_key
                .smart_mul_parallelized(&mut u, &mut a_ge_b)
                .clone(),
        );
        let sv_c2 = server_key.smart_add_parallelized(
            &mut server_key
                .smart_mul_parallelized(&mut sv, &mut a_lt_b)
                .clone(),
            &mut server_key
                .smart_mul_parallelized(&mut v, &mut a_ge_b)
                .clone(),
        );
        let a_c2 = server_key
            .scalar_right_shift_parallelized(&server_key.sub_parallelized(&sa_c2, &sb_c2), 1);

        let u_s_v = sub_mod::<NB, _>(&su_c2, &sv_c2, p, server_key);
        let u_s_v_d2 = mul_mod::<NB, _>(
            &u_s_v,
            &server_key.create_trivial_radix(126u64, NB), // TODO: fix this
            p,
            server_key,
        );
        let u_c2 = mersenne::mod_mersenne::<NB, _>(&u_s_v_d2, p, server_key);
        let b_c2 = sb_c2;
        let v_c2 = sv_c2;

        // consolidate c1 and c2

        a = server_key.smart_add_parallelized(
            &mut server_key
                .smart_mul_parallelized(&mut a_c1.clone(), &mut a_is_even.clone())
                .clone(),
            &mut server_key
                .smart_mul_parallelized(&mut a_c2.clone(), &mut a_is_odd.clone())
                .clone(),
        );
        b = server_key.smart_add_parallelized(
            &mut server_key
                .smart_mul_parallelized(&mut b_c1.clone(), &mut a_is_even.clone())
                .clone(),
            &mut server_key
                .smart_mul_parallelized(&mut b_c2.clone(), &mut a_is_odd.clone())
                .clone(),
        );
        u = server_key.smart_add_parallelized(
            &mut server_key
                .smart_mul_parallelized(&mut u_c1.clone(), &mut a_is_even.clone())
                .clone(),
            &mut server_key
                .smart_mul_parallelized(&mut u_c2.clone(), &mut a_is_odd.clone())
                .clone(),
        );
        v = server_key.smart_add_parallelized(
            &mut server_key
                .smart_mul_parallelized(&mut v_c1.clone(), &mut a_is_even.clone())
                .clone(),
            &mut server_key
                .smart_mul_parallelized(&mut v_c2.clone(), &mut a_is_odd.clone())
                .clone(),
        );
        let mut done = server_key.smart_scalar_eq_parallelized(&mut a.clone(), 0);
        let mut never_done = server_key
            .smart_sub_parallelized(&mut server_key.create_trivial_radix(1, 1), &mut was_done);
        let mut done_now = server_key.smart_bitand_parallelized(&mut done, &mut never_done);
        server_key.smart_bitor_assign_parallelized(&mut was_done, &mut done);
        // inv = inv + done_now * t1
        let mut update = server_key.smart_mul_parallelized(&mut done_now, &mut v.clone());
        server_key.smart_add_assign_parallelized(&mut result, &mut update);
        read_client_key(|server_key| {
            println!(
                "a = {}\nb = {}\nu = {}\nv = {}\nresult {}\n-----",
                format(server_key.decrypt_radix::<P>(&a)),
                format(server_key.decrypt_radix::<P>(&b)),
                format(server_key.decrypt_radix::<P>(&u)),
                format(server_key.decrypt_radix::<P>(&v)),
                format(server_key.decrypt_radix::<P>(&result))
            )
        });
    }
    result.clone()
}

/// a^-1 mod p where a*a^-1 = 1 mod p
#[inline]
pub fn inverse_mod<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    inverse_mod_trim::<NB, _>(a, p, server_key)
}

/// a^-1 mod p where a*a^-1 = 1 mod p
pub fn inverse_mod_trim<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    #[cfg(feature = "low_level_timing")]
    let start_ops = Instant::now();
    #[cfg(feature = "low_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);
    #[cfg(feature = "low_level_timing")]
    println!("inverse_mod: start -- ref {}", task_ref);

    let padded_nb = NB + 1;
    // implement extended euclidean algorithm with trim bit
    // assume a < p. (no check)
    let a = server_key.extend_radix_with_trivial_zero_blocks_msb(&a.clone(), 1);
    let mut r0 = server_key.create_trivial_radix(p, padded_nb);
    let mut r1 = a;
    let mut was_done = server_key.create_trivial_radix(0, 1);
    let mut t0 = server_key.create_trivial_radix(0, padded_nb);
    let mut t1 = server_key.create_trivial_radix(1, padded_nb);
    let mut inv = server_key.create_trivial_radix(0, padded_nb);
    let mut trim = 0;
    // euclidean algorithm
    // NB/2 best case and NB worst case
    let loop_end = <P as Numeric>::BITS + 1;
    for i in 0..loop_end {
        #[cfg(feature = "low_level_timing")]
        let bit_start = Instant::now();
        println!("inverse_mod: bit {}", i);

        // q, r = r0 / r1
        let (mut q, r) = server_key.smart_div_rem_parallelized(&mut r0.clone(), &mut r1.clone());

        // rayon::join(
        //     || server_key.full_propagate_parallelized(&mut q),
        //     || server_key.full_propagate_parallelized(&mut r),
        // );
        server_key.extend_radix_with_trivial_zero_blocks_msb_assign(&mut q, trim);
        let full_r = server_key.extend_radix_with_trivial_zero_blocks_msb(&r, trim);

        let tmp = t1.clone();

        let mut qt1 = server_key.smart_mul_parallelized(&mut t1, &mut q);
        // t1 = t0 - q * t1
        t1 = server_key.smart_sub_parallelized(&mut t0, &mut qt1);
        t0 = tmp;
        // is_done = r =? 0
        // never_done = 1 - is_done
        // was_done = was_done | is_done
        // done_now = is_done & never_done
        let mut done = server_key.smart_scalar_eq_parallelized(&mut full_r.clone(), 0);
        let len = done.blocks().len();
        server_key.trim_radix_blocks_msb_assign(&mut done, len - 1);
        let mut never_done = server_key
            .smart_sub_parallelized(&mut server_key.create_trivial_radix(1, 1), &mut was_done);
        let done_now = server_key.smart_bitand_parallelized(&mut done, &mut never_done);
        server_key.smart_bitor_assign_parallelized(&mut was_done, &mut done);

        let len = done_now.blocks().len();
        let mut not_done_now = server_key.smart_sub_parallelized(
            &mut server_key.create_trivial_radix(0, padded_nb),
            &mut server_key.extend_radix_with_trivial_zero_blocks_msb(&done_now, padded_nb - len),
        );

        let mut update = server_key.smart_bitand_parallelized(&mut t0, &mut not_done_now);
        //let mut update = server_key.smart_mul_parallelized(&mut t0, &mut done_now);
        server_key.smart_add_assign_parallelized(&mut inv, &mut update);

        // update values
        if (i % 2 == 0) & (i != 0) {
            r0 = server_key.trim_radix_blocks_msb(&r1.clone(), 1);
            r1 = server_key.trim_radix_blocks_msb(&r.clone(), 1);
            trim += 1;
        } else {
            r0 = r1.clone();
            r1 = r.clone();
        }

        #[cfg(feature = "low_level_timing")]
        {
            if (i == 0) | (i == loop_end - 1) {
                println!(
                    "----Inverse mod bit {i} done in {:.2}s -- ref {}",
                    bit_start.elapsed().as_secs_f64(),
                    task_ref
                );
            }
        }
    }

    // final result mod p
    // inverse can be **negative**. so we need to add p to make it positive

    server_key.smart_scalar_add_assign_parallelized(&mut inv, p);

    let mut is_gt = server_key.smart_scalar_ge_parallelized(&mut inv, p);
    server_key.trim_radix_blocks_msb_assign(&mut is_gt, padded_nb - 1);

    let mut to_sub = server_key.smart_mul_parallelized(
        &mut server_key.create_trivial_radix(p, padded_nb),
        &mut is_gt,
    );
    server_key.smart_sub_assign_parallelized(&mut inv, &mut to_sub);
    server_key.full_propagate_parallelized(&mut inv);

    #[cfg(feature = "low_level_timing")]
    println!(
        "Inverse mod done in {:.2}s -- ref {}",
        start_ops.elapsed().as_secs_f64(),
        task_ref
    );

    inv
}

/// a^-1 mod p where a*a^-1 = 1 mod p
pub fn inverse_mod_without_trim<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let padded_nb = NB + 1;
    // implement extended euclidean algorithm
    // assume a < p. (no check)
    let a = server_key.extend_radix_with_trivial_zero_blocks_msb(&a.clone(), 1);
    let mut r0 = server_key.create_trivial_radix(p, padded_nb);
    let mut r1 = a;
    let mut was_done = server_key.create_trivial_radix(0, 1);
    let mut t0 = server_key.create_trivial_radix(0, padded_nb);
    let mut t1 = server_key.create_trivial_radix(1, padded_nb);
    let mut inv = server_key.create_trivial_radix(0, padded_nb);

    // euclidean algorithm
    // NB/2 best case and NB worst case
    for _i in 0..<P as Numeric>::BITS {
        let _now = Instant::now();
        // q, r = r0 / r1
        let (q, r) = server_key.smart_div_rem_parallelized(&mut r0.clone(), &mut r1.clone());
        // rayon::join(
        //     || server_key.full_propagate_parallelized(&mut q),
        //     || server_key.full_propagate_parallelized(&mut r),
        // );
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
    }

    // final result mod p
    // inverse can be **negative**. so we need to add p to make it positive
    server_key.smart_scalar_add_assign_parallelized(&mut inv, p);
    let mut is_gt = server_key.smart_scalar_ge_parallelized(&mut inv, p);
    server_key.trim_radix_blocks_msb_assign(&mut is_gt, NB - 1);
    let mut to_sub =
        server_key.smart_mul_parallelized(&mut server_key.create_trivial_radix(p, NB), &mut is_gt);
    server_key.smart_sub_assign_parallelized(&mut inv, &mut to_sub);
    // server_key.full_propagate_parallelized(&mut inv);
    inv
}

/// a + b mod p
pub fn add_mod<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    // assume large p and a,b < p
    let start_ops = Instant::now();
    #[cfg(feature = "low_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);

    #[cfg(feature = "low_level_timing")]
    println!("Add mod start -- ref {}", task_ref);

    let mut a_expanded = server_key.extend_radix_with_trivial_zero_blocks_msb(a, 1);
    server_key.smart_add_assign_parallelized(&mut a_expanded, &mut b.clone());
    let res = modulo_fast::<NB, _>(&a_expanded, p, server_key);

    ProtocolStats::add_time(ProtocolLowOps::AddMod, start_ops.elapsed().as_secs_f32());
    #[cfg(feature = "low_level_timing")]
    println!(
        "Add mod done in {:.2}s -- ref {}",
        start_ops.elapsed().as_secs_f64(),
        task_ref
    );

    res
}

/// a - b mod p
pub fn sub_mod<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let start_ops = Instant::now();
    #[cfg(feature = "low_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);
    #[cfg(feature = "low_level_timing")]
    println!("Sub mod start -- ref {}", task_ref);

    let is_gt = server_key.smart_gt_parallelized(&mut b.clone(), &mut a.clone());
    let mut to_add = selector_zero_constant::<NB, _>(p, &is_gt, server_key);
    let mut a_expanded = server_key.extend_radix_with_trivial_zero_blocks_msb(a, 1);
    server_key.smart_add_assign_parallelized(&mut a_expanded, &mut to_add);
    server_key.sub_assign_parallelized(&mut a_expanded, b);
    //server_key.full_propagate_parallelized(&mut a_expanded);
    server_key.trim_radix_blocks_msb_assign(&mut a_expanded, 1);

    ProtocolStats::add_time(ProtocolLowOps::SubMod, start_ops.elapsed().as_secs_f32());
    #[cfg(feature = "low_level_timing")]
    println!(
        "Sub mod done in {:.2}s -- ref {}",
        start_ops.elapsed().as_secs_f64(),
        task_ref
    );
    a_expanded
}

/// a * b mod p
pub fn mul_mod_bitwise<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    // assume large p and a,b < p
    #[cfg(feature = "low_level_timing")]
    let start_ops = Instant::now();
    #[cfg(feature = "low_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);

    #[cfg(feature = "low_level_timing")]
    println!("Mul mod bitwise start -- ref {}", task_ref);

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

    let loop_end = <P as Numeric>::BITS;
    for _i in 0..loop_end {
        #[cfg(feature = "low_level_timing")]
        let bit_start = Instant::now();

        ((b_next_tmp, a_tmp), (bit, (res, to_add_later))) = rayon::join(
            || {
                rayon::join(
                    || server_key.scalar_right_shift_parallelized(&b_next_tmp, 1),
                    || double_mod::<NB, _>(&a_tmp, p, server_key),
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
        #[cfg(feature = "low_level_timing")]
        {
            if (_i == 0) | (_i == loop_end - 1) {
                println!(
                    "----Mul mod bitwise {_i} done in {:.2}s -- ref {}",
                    bit_start.elapsed().as_secs_f64(),
                    task_ref
                );
            }
        }
    }

    let result = add_mod::<NB, _>(&res, &to_add_later, p, server_key);

    #[cfg(feature = "low_level_timing")]
    println!(
        "Mul mod bitwise done in {:.2}s -- ref {}",
        start_ops.elapsed().as_secs_f64(),
        task_ref
    );

    result
}

/// a * b mod p
pub fn mul_mod_div_rem<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    // assume large p and a,b < p
    #[cfg(feature = "low_level_timing")]
    let start_ops = Instant::now();
    #[cfg(feature = "low_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);
    #[cfg(feature = "low_level_timing")]
    println!("Mul mod div rem start -- ref {}", task_ref);

    let mut a_expanded = server_key.extend_radix_with_trivial_zero_blocks_msb(a, NB);
    server_key.smart_mul_assign_parallelized(&mut a_expanded, &mut b.clone());
    // server_key.full_propagate_parallelized(&mut a_expanded);
    let (_q, mut r) = server_key.smart_div_rem_parallelized(
        &mut a_expanded,
        &mut server_key.create_trivial_radix(p, NB * 2),
    );
    // server_key.full_propagate_parallelized(&mut r);
    server_key.trim_radix_blocks_msb_assign(&mut r, NB);
    #[cfg(feature = "low_level_timing")]
    println!(
        "Mul mod div rem done in {:.2}s -- ref {}",
        start_ops.elapsed().as_secs_f64(),
        task_ref
    );
    r
}

/// a * b mod p
pub fn mul_mod<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let start_ops = Instant::now();
    let res = mul_mod_mersenne::<NB, _>(a, b, p, server_key);
    ProtocolStats::add_time(ProtocolLowOps::MulMod, start_ops.elapsed().as_secs_f32());
    res
}

/// a * b mod p where b is a constant
/// slower than 12 `add_mod`
pub fn mul_mod_constant<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    b: P,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    // assume large p and a,b < p
    #[cfg(feature = "low_level_timing")]
    let start_ops = Instant::now();
    #[cfg(feature = "low_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);
    #[cfg(feature = "low_level_timing")]
    println!("Mul mod constant start -- ref {}", task_ref);

    let mut a_expanded = server_key.extend_radix_with_trivial_zero_blocks_msb(a, NB);
    server_key.smart_scalar_mul_assign_parallelized(&mut a_expanded, b);
    let res = mod_mersenne::<NB, _>(&a_expanded, p, server_key);
    #[cfg(feature = "low_level_timing")]
    println!(
        "Mul mod constant done in {:.2}s -- ref {}",
        start_ops.elapsed().as_secs_f64(),
        task_ref
    );
    res
}

/// a^2 mod p
#[inline(always)]
pub fn square_mod<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let start_ops = Instant::now();
    let res = mul_mod_mersenne::<NB, _>(a, a, p, server_key);
    ProtocolStats::add_time(ProtocolLowOps::SquareMod, start_ops.elapsed().as_secs_f32());
    res
}

/// a*2 mod p
#[inline(always)]
pub fn double_mod<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let start_ops = Instant::now();
    #[cfg(feature = "low_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);
    #[cfg(feature = "low_level_timing")]
    println!("Double mod start -- ref {}", task_ref);

    let mut a_expanded = server_key.extend_radix_with_trivial_zero_blocks_msb(a, 1);
    server_key.scalar_left_shift_assign_parallelized(&mut a_expanded, 1);
    let res = modulo_fast::<NB, _>(&a_expanded, p, server_key);

    ProtocolStats::add_time(ProtocolLowOps::DoubleMod, start_ops.elapsed().as_secs_f32());
    #[cfg(feature = "low_level_timing")]
    println!(
        "Double mod done in {:.2}s -- ref {}",
        start_ops.elapsed().as_secs_f64(),
        task_ref
    );
    res
}

/// a^b mod p
pub fn pow_mod<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    #[cfg(feature = "low_level_timing")]
    let start_ops = Instant::now();
    #[cfg(feature = "low_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);
    #[cfg(feature = "low_level_timing")]
    println!("Pow mod start -- ref {}", task_ref);

    let mut res = server_key.create_trivial_radix(1, NB);
    let mut base = a.clone();
    let mut exponent = b.clone();
    let loop_end = <P as Numeric>::BITS;
    for _i in 0..loop_end {
        #[cfg(feature = "low_level_timing")]
        let bit_start = Instant::now();

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
        #[cfg(feature = "low_level_timing")]
        {
            if (_i == 0) | (_i == loop_end - 1) {
                println!(
                    "----pow mod bit {_i} done in {:.2}s -- ref {}",
                    bit_start.elapsed().as_secs_f32(),
                    task_ref,
                );
            }
        }
    }
    #[cfg(feature = "low_level_timing")]
    println!(
        "Pow mod done in {:.2}s -- ref {}",
        start_ops.elapsed().as_secs_f64(),
        task_ref
    );
    res
}

pub fn inverse_mod_pow<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    #[cfg(feature = "low_level_timing")]
    let start_ops = Instant::now();
    #[cfg(feature = "low_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);
    #[cfg(feature = "low_level_timing")]
    println!("Inverse mod start -- ref {}", task_ref);

    let mut res = server_key.create_trivial_radix(1, NB);
    let mut base = a.clone();
    let mut exponent = p - P::cast_from(2u64);
    let loop_end = <P as Numeric>::BITS;
    for _i in 0..loop_end {
        #[cfg(feature = "low_level_timing")]
        let bit_start = Instant::now();
        let b = base.clone();
        rayon::join(
            || {
                if exponent & P::ONE == P::ONE {
                    res = mul_mod::<NB, _>(&res, &b, p, server_key);
                }
                exponent >>= 1;
            },
            || base = square_mod::<NB, _>(&base, p, server_key),
        );
        #[cfg(feature = "low_level_timing")]
        {
            if (_i == 0) | (_i == loop_end - 1) {
                println!(
                    "----inverse mod bit {_i} done in {:.2}s -- ref {}",
                    bit_start.elapsed().as_secs_f32(),
                    task_ref,
                );
            }
        }
    }
    #[cfg(feature = "low_level_timing")]
    println!(
        "Inverse mod done in {:.2}s -- ref {}",
        start_ops.elapsed().as_secs_f64(),
        task_ref
    );
    res
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use tfhe::{
        integer::{keycache::IntegerKeyCache, U256},
        shortint::prelude::PARAM_MESSAGE_2_CARRY_2,
    };

    use crate::{
        numeral::Numeral,
        ops::{
            add_mod, double_mod, inverse_mod,
            mersenne::mod_mersenne,
            mul_mod, mul_mod_constant, multi_add_mod,
            native::{inverse_mod_native, mul_mod_native, sub_mod_native},
            sub_mod,
        },
        CLIENT_KEY,
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
        println!("8a using addition - {:.2}s", timer.elapsed().as_secs());

        let timer = Instant::now();
        let _enc_8a_mul = mul_mod_constant::<NUM_BLOCK, _>(&enc_a, 8, p, &server_key);
        println!(
            "8a using multiplication - {:.2}s",
            timer.elapsed().as_secs()
        );
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
        let p: u128 = 251;
        let a: u128 = 248;
        let b: u128 = 249;
        let c: u128 = sub_mod_native(a, b, p);

        let enc_c = sub_mod::<NUM_BLOCK, _>(
            &client_key.encrypt_radix(a, NUM_BLOCK),
            &client_key.encrypt_radix(b, NUM_BLOCK),
            p,
            &server_key,
        );
        assert_eq!(c as u8, client_key.decrypt_radix::<u8>(&enc_c));

        let d = sub_mod_native(c, b, p);
        let enc_d = sub_mod::<NUM_BLOCK, _>(
            &enc_c,
            &client_key.encrypt_radix(b, NUM_BLOCK),
            p,
            &server_key,
        );
        assert_eq!(d as u8, client_key.decrypt_radix::<u8>(&enc_d));

        let e = sub_mod_native(c, d, p);
        let enc_e = sub_mod::<NUM_BLOCK, _>(&enc_c, &enc_d, p, &server_key);
        assert_eq!(e as u8, client_key.decrypt_radix::<u8>(&enc_e));
    }

    #[test]
    fn correct_mul_mod() {
        let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);
        *CLIENT_KEY.write().unwrap() = Some(client_key.clone());
        const NUM_BLOCK: usize = 4;
        let p: u128 = 251;
        let a: u128 = 249;
        let b: u128 = 248;
        let c: u128 = mul_mod_native(a, b, p);

        let enc_c = mul_mod::<NUM_BLOCK, _>(
            &client_key.encrypt_radix(a, NUM_BLOCK),
            &client_key.encrypt_radix(b, NUM_BLOCK),
            p,
            &server_key,
        );
        assert_eq!(c as u8, client_key.decrypt_radix::<u8>(&enc_c));

        let d = mul_mod_native(c, b, p);
        let enc_d = mul_mod::<NUM_BLOCK, _>(
            &enc_c,
            &client_key.encrypt_radix(b, NUM_BLOCK),
            p,
            &server_key,
        );
        assert_eq!(d as u8, client_key.decrypt_radix::<u8>(&enc_d));

        let e = mul_mod_native(c, d, p);
        let enc_e = mul_mod::<NUM_BLOCK, _>(&enc_c, &enc_d, p, &server_key);
        assert_eq!(e as u8, client_key.decrypt_radix::<u8>(&enc_e));

        let f = mul_mod_native(e, e, p);
        let enc_f = mul_mod::<NUM_BLOCK, _>(&enc_e, &enc_e, p, &server_key);
        assert_eq!(f as u8, client_key.decrypt_radix::<u8>(&enc_f));

        let g = mul_mod_native(f, f, p);
        let enc_g = mul_mod::<NUM_BLOCK, _>(&enc_f, &enc_f, p, &server_key);
        assert_eq!(g as u8, client_key.decrypt_radix::<u8>(&enc_g));
    }

    #[test]
    fn correct_multi_add_mod() {
        let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);
        const NUM_BLOCK: usize = 4;
        let p: u8 = 251;

        let multi_add_mod_naive =
            |a: &[u128]| a.iter().copied().fold(0, |a, b| (a + b) % p as u128);

        let a: [u128; 4] = [248, 249, 250, 251];
        let b: [u128; 4] = [251, 250, 249, 248];
        let c: u128 = multi_add_mod_naive(&a);
        let enc_c = multi_add_mod::<NUM_BLOCK, _>(
            &a.iter()
                .map(|a| client_key.encrypt_radix(*a, NUM_BLOCK))
                .collect::<Vec<_>>(),
            p,
            &server_key,
        );
        assert_eq!(c as u8, client_key.decrypt_radix::<u8>(&enc_c));

        let d = multi_add_mod_naive(&b);
        let enc_d = multi_add_mod::<NUM_BLOCK, _>(
            &b.iter()
                .map(|a| client_key.encrypt_radix(*a, NUM_BLOCK))
                .collect::<Vec<_>>(),
            p,
            &server_key,
        );
        assert_eq!(d as u8, client_key.decrypt_radix::<u8>(&enc_d));
    }

    #[test]
    fn correct_inverse_mod() {
        let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);
        const NUM_BLOCK: usize = 4;
        let p: u128 = 157;
        let a: u128 = 8;
        let b: u128 = 6;
        let e: u128 = 45;
        let f: u128 = 123;
        let h: u128 = 127;
        let i: u128 = 156;

        let c = inverse_mod_native(a, p);
        let enc_c =
            inverse_mod::<NUM_BLOCK, _>(&client_key.encrypt_radix(a, NUM_BLOCK), p, &server_key);
        assert_eq!(c as u8, client_key.decrypt_radix::<u8>(&enc_c));

        let d = inverse_mod_native(b, p);
        let enc_d =
            inverse_mod::<NUM_BLOCK, _>(&client_key.encrypt_radix(b, NUM_BLOCK), p, &server_key);
        assert_eq!(d as u8, client_key.decrypt_radix::<u8>(&enc_d));

        let m = inverse_mod_native(e, p);
        let enc_m =
            inverse_mod::<NUM_BLOCK, _>(&client_key.encrypt_radix(e, NUM_BLOCK), p, &server_key);
        assert_eq!(m as u8, client_key.decrypt_radix::<u8>(&enc_m));

        let g = inverse_mod_native(f, p);
        let enc_g =
            inverse_mod::<NUM_BLOCK, _>(&client_key.encrypt_radix(f, NUM_BLOCK), p, &server_key);
        assert_eq!(g as u8, client_key.decrypt_radix::<u8>(&enc_g));

        let j = inverse_mod_native(h, p);
        let enc_j =
            inverse_mod::<NUM_BLOCK, _>(&client_key.encrypt_radix(h, NUM_BLOCK), p, &server_key);
        assert_eq!(j as u8, client_key.decrypt_radix::<u8>(&enc_j));

        let k = inverse_mod_native(i, p);
        let enc_k =
            inverse_mod::<NUM_BLOCK, _>(&client_key.encrypt_radix(i, NUM_BLOCK), p, &server_key);
        assert_eq!(k as u8, client_key.decrypt_radix::<u8>(&enc_k));

        let l = inverse_mod_native(f, p);
        let enc_l =
            inverse_mod::<NUM_BLOCK, _>(&client_key.encrypt_radix(f, NUM_BLOCK), p, &server_key);
        assert_eq!(l as u8, client_key.decrypt_radix::<u8>(&enc_l));
    }

    #[test]
    fn parallel_mod_reduc_mul() {
        let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);

        const NUM_BLOCK: usize = 4;
        type Integer = u8;
        let p: Integer = 251;
        let x1: Integer = 8;
        let y1: Integer = 45;

        let ct_x1 = client_key.encrypt_radix(x1, NUM_BLOCK);
        let ct_y1 = client_key.encrypt_radix(y1, NUM_BLOCK);

        let now = Instant::now();
        let _ = mul_mod::<NUM_BLOCK, _>(&ct_x1, &ct_y1, p, &server_key);
        println!("mul_mod: {:.2}s", now.elapsed().as_secs_f32());

        let now = Instant::now();
        let _ = rayon::join(
            || {
                let mut expanded =
                    server_key.extend_radix_with_trivial_zero_blocks_msb(&ct_x1, NUM_BLOCK);
                server_key.smart_mul_assign_parallelized(&mut expanded, &mut ct_y1.clone());
                expanded
            },
            || mod_mersenne::<NUM_BLOCK, _>(&ct_x1, p, &server_key),
        );
        println!("parallel_mod_mersenne: {:.2}s", now.elapsed().as_secs_f32());
    }

    #[test]
    fn correct_sub_wrap() {
        let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);

        const NUM_BLOCK: usize = 8;

        let x1 = 1u8;
        let mut ct_x1 = client_key.encrypt_radix(x1, 1);

        let res = server_key.smart_sub_parallelized(
            &mut server_key.create_trivial_radix(0, NUM_BLOCK + 1),
            &mut ct_x1,
        );

        println!("res: {}", u8::decrypt_bigint(&res, &client_key));
    }
}
