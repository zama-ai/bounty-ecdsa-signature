use std::time::Instant;

use logging_timer::{stime, time, timer, Level};
use rand::Rng;
use rayon::prelude::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use tfhe::{
    core_crypto::prelude::{Numeric, UnsignedInteger},
    integer::{
        block_decomposition::{DecomposableInto, RecomposableFrom},
        server_key::MiniUnsignedInteger,
        IntegerCiphertext, RadixCiphertext, ServerKey,
    },
};

use crate::{
    helper::{format, read_client_key},
    numeral::Numeral,
    ops::mersenne::mod_mersenne,
    stats::{ProtocolLowOps, ProtocolStats},
};

use self::{mersenne::mul_mod_mersenne, native::inverse_mod_native, primitive::parallel_fn};

pub mod group_homogenous;
pub mod group_jacobian;
pub mod mersenne;
pub mod native;
pub mod primitive;
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
    #[cfg(not(feature = "inw_selector"))]
    {
        let mut selector = server_key.trim_radix_blocks_msb(&selector, len - 1);
        server_key.smart_mul_assign_parallelized(&mut res, &mut selector);
    }
    #[cfg(feature = "inw_selector")]
    {
        let a_len = a.blocks().len();
        let mut not_selector = server_key.smart_sub_parallelized(
            &mut server_key.create_trivial_radix(0, a_len),
            &mut server_key.extend_radix_with_trivial_zero_blocks_msb(selector, a_len - len),
        );
        server_key.smart_bitand_assign_parallelized(&mut res, &mut not_selector);
    }
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
            &mut server_key.extend_radix_with_trivial_zero_blocks_msb(selector, NB - len),
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
    not_selector: &RadixCiphertext,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let selector_len = selector.blocks().len();
    let not_selector_len = not_selector.blocks().len();
    let mut selector = server_key.trim_radix_blocks_msb(selector, selector_len - 1);
    let mut not_selector = server_key.trim_radix_blocks_msb(not_selector, not_selector_len - 1);
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
    // add all elements in a together
    let le = a.len();
    // find bit length of le
    let log_le = UnsignedInteger::ceil_ilog2(le) as usize;
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
    sum
}

/// turn x mod a to x mod b
/// only if a > b and a < 2b
#[time("trace", "Modulus Reduction")]
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
pub fn modulo_div_rem<const NB: usize, P: Numeral>(
    x: &RadixCiphertext,
    b: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let (_q, r) = server_key.smart_div_rem_parallelized(
        &mut x.clone(),
        &mut server_key.create_trivial_radix(b, x.blocks().len()),
    );
    r
}

pub fn inverse_mod_binary_gcd<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    // implement binary gcd algorithm
    // assume a < p. (no check)
    let div_two = inverse_mod_native(P::TWO, p);

    let mul_mod_div_two = |a: &RadixCiphertext| {
        let mut res = server_key.smart_scalar_mul_parallelized(
            &mut server_key.extend_radix_with_trivial_zero_blocks_msb(a, NB),
            div_two,
        );
        mod_mersenne::<NB, _>(&mut res, p, server_key)
    };

    let mut a = a.clone();
    let mut b = server_key.create_trivial_radix(p, NB);
    let mut u: RadixCiphertext = server_key.create_trivial_radix(1, NB);
    let mut v: RadixCiphertext = server_key.create_trivial_radix(0, NB);

    read_client_key(|server_key| {
        println!(
            "a = {}, b = {}, u = {}, v = {} -----",
            format(server_key.decrypt_radix::<P>(&a)),
            format(server_key.decrypt_radix::<P>(&b)),
            format(server_key.decrypt_radix::<P>(&u)),
            format(server_key.decrypt_radix::<P>(&v)),
        )
    });

    // if condition return a,b else b,a
    let select_two = |condition: &RadixCiphertext, a: &RadixCiphertext, b: &RadixCiphertext| {
        let (mut added, mut sel_a) = rayon::join(
            || {
                server_key.smart_add_parallelized(
                    &mut server_key.extend_radix_with_trivial_zero_blocks_msb(a, 1),
                    &mut b.clone(),
                )
            },
            || {
                server_key.if_then_else_parallelized(
                    &mut condition.clone(),
                    &mut a.clone(),
                    &mut b.clone(),
                )
            },
        );
        let sel_b = server_key.smart_sub_parallelized(&mut added, &mut sel_a);
        (sel_a, server_key.trim_radix_blocks_msb(&sel_b, 1))
    };

    for _i in 0..<P as Numeric>::BITS * 2 - 1 {
        let _tmr = timer!(Level::Trace; "Inverse Mod Binary GCD", "Bit {}", _i);

        let (mut a_mod_two, mut a_lt_b) = rayon::join(
            || server_key.smart_scalar_bitand_parallelized(&mut a.clone(), 1),
            || server_key.smart_lt_parallelized(&mut a.clone(), &mut b),
        );
        server_key.trim_radix_blocks_msb_assign(&mut a_mod_two, NB - 1);
        server_key.trim_radix_blocks_msb_assign(&mut a_lt_b, NB - 1);
        let a_mod_two_and_lt_b = server_key.smart_bitand_parallelized(&mut a_mod_two, &mut a_lt_b);

        ((a, b), (u, v)) = rayon::join(
            || select_two(&a_mod_two_and_lt_b, &b, &a),
            || select_two(&a_mod_two_and_lt_b, &v, &u),
        );

        (a, u) = rayon::join(
            || {
                let mult = server_key.smart_mul_parallelized(&mut b, &mut a_mod_two.clone());
                sub_mod::<NB, _>(&a, &mult, p, server_key)
                //server_key.sub_parallelized(&a, &mult)
            },
            || {
                let mult = server_key.smart_mul_parallelized(&mut v, &mut a_mod_two.clone());
                sub_mod::<NB, _>(&u, &mult, p, server_key)
                //server_key.sub_parallelized(&u, &mult)
            },
        );

        (a, u) = rayon::join(
            || server_key.scalar_right_shift_parallelized(&mut a, 1),
            || mul_mod_div_two(&u),
        );

        read_client_key(|server_key| {
            println!(
                "a = {}, b = {}, u = {}, v = {} -----",
                format(server_key.decrypt_radix::<P>(&a)),
                format(server_key.decrypt_radix::<P>(&b)),
                format(server_key.decrypt_radix::<P>(&u)),
                format(server_key.decrypt_radix::<P>(&v)),
            )
        });
    }

    v
}

/// a^-1 mod p where a*a^-1 = 1 mod p
#[inline]
#[time("debug", "Inverse Mod")]
pub fn inverse_mod<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    #[cfg(feature = "binary_gcd")]
    return inverse_mod_binary_gcd::<NB, _>(a, p, server_key);
    #[cfg(not(feature = "binary_gcd"))]
    return inverse_mod_trim::<NB, _>(a, p, server_key);
}

#[inline]
pub fn inverse_mods<const NB: usize, P: Numeral>(
    a: &[RadixCiphertext],
    p: P,
    server_key: &ServerKey,
) -> Vec<RadixCiphertext> {
    let product = parallel_fn(a, |a, b| mul_mod::<NB, _>(a, b, p, server_key));
    let inversed = inverse_mod::<NB, _>(&product, p, server_key);
    let mut result = vec![server_key.create_trivial_radix(0, NB); a.len()];

    (0..a.len())
        .into_par_iter()
        .map(|i| {
            let mut coef = a
                .iter()
                .cloned()
                .enumerate()
                .filter(|(j, _)| i != *j)
                .map(|e| e.1)
                .collect::<Vec<_>>();
            coef.push(inversed.clone());
            parallel_fn(&coef, |a, b| mul_mod::<NB, _>(a, b, p, server_key))
        })
        .collect_into_vec(&mut result);

    result
}

/// a^-1 mod p where a*a^-1 = 1 mod p
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
        let (mut q, r) = server_key.smart_div_rem_parallelized(&mut r0.clone(), &mut r1.clone());

        //read_client_key(|ck| {
        //println!("q: {}", P::decrypt(&q, &ck).format());
        //println!("r: {}", P::decrypt(&r, &ck).format());
        //println!("blocks r: {}, q: {}", r.blocks().len(), q.blocks().len());
        //});

        //rayon::join(
        //|| server_key.full_propagate_parallelized(&mut q),
        //|| server_key.full_propagate_parallelized(&mut r),
        //);
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

        //let len = done_now.blocks().len();
        //let mut not_done_now = server_key.smart_sub_parallelized(
        //&mut server_key.create_trivial_radix(0, padded_nb),
        //&mut server_key.extend_radix_with_trivial_zero_blocks_msb(&done_now, padded_nb - len),
        //);

        //let mut update = server_key.smart_bitand_parallelized(&mut t0, &mut not_done_now);
        //let mut update = server_key.smart_mul_parallelized(&mut t0, &mut done_now);
        let mut update = selector_zero(&t0, &done_now, server_key);
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
    let mut r0: RadixCiphertext = server_key.create_trivial_radix(p, padded_nb);
    let mut r1 = a;
    let mut was_done = server_key.create_trivial_radix(0, 1);
    let mut t0: RadixCiphertext = server_key.create_trivial_radix(0, padded_nb);
    let mut t1: RadixCiphertext = server_key.create_trivial_radix(1, padded_nb);
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
#[time("debug", "Add Mod")]
pub fn add_mod<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    // assume large p and a,b < p
    let start_ops = Instant::now();

    let mut a_expanded = server_key.extend_radix_with_trivial_zero_blocks_msb(a, 1);
    server_key.smart_add_assign_parallelized(&mut a_expanded, &mut b.clone());
    let res = modulo_fast::<NB, _>(&a_expanded, p, server_key);

    ProtocolStats::add_time(ProtocolLowOps::AddMod, start_ops.elapsed().as_secs_f32());

    res
}

/// a - b mod p
#[time("debug", "Sub Mod")]
pub fn sub_mod<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let start_ops = Instant::now();

    let is_gt = server_key.smart_gt_parallelized(&mut b.clone(), &mut a.clone());
    let mut to_add = selector_zero_constant::<NB, _>(p, &is_gt, server_key);
    let mut a_expanded = server_key.extend_radix_with_trivial_zero_blocks_msb(a, 1);
    server_key.smart_add_assign_parallelized(&mut a_expanded, &mut to_add);
    server_key.sub_assign_parallelized(&mut a_expanded, b);
    //server_key.full_propagate_parallelized(&mut a_expanded);
    server_key.trim_radix_blocks_msb_assign(&mut a_expanded, 1);

    ProtocolStats::add_time(ProtocolLowOps::SubMod, start_ops.elapsed().as_secs_f32());
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
    let mut res: RadixCiphertext = server_key.create_trivial_radix(0u64, NB);
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
    }

    add_mod::<NB, _>(&res, &to_add_later, p, server_key)
}

/// a * b mod p
pub fn mul_mod_div_rem<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    // assume large p and a,b < p
    let mut a_expanded = server_key.extend_radix_with_trivial_zero_blocks_msb(a, NB);
    server_key.smart_mul_assign_parallelized(&mut a_expanded, &mut b.clone());
    // server_key.full_propagate_parallelized(&mut a_expanded);
    let (_q, mut r) = server_key.smart_div_rem_parallelized(
        &mut a_expanded,
        &mut server_key.create_trivial_radix(p, NB * 2),
    );
    // server_key.full_propagate_parallelized(&mut r);
    server_key.trim_radix_blocks_msb_assign(&mut r, NB);
    r
}

/// a * b mod p
#[time("debug", "Mul Mod")]
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
    let mut a_expanded = server_key.extend_radix_with_trivial_zero_blocks_msb(a, NB);
    server_key.smart_scalar_mul_assign_parallelized(&mut a_expanded, b);

    mod_mersenne::<NB, _>(&a_expanded, p, server_key)
}

/// a^2 mod p
#[inline(always)]
#[time("debug", "Square mod")]
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
#[time("debug", "Double mod")]
pub fn double_mod<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let start_ops = Instant::now();
    let mut a_expanded = server_key.extend_radix_with_trivial_zero_blocks_msb(a, 1);
    server_key.scalar_left_shift_assign_parallelized(&mut a_expanded, 1);
    let res = modulo_fast::<NB, _>(&a_expanded, p, server_key);
    ProtocolStats::add_time(ProtocolLowOps::DoubleMod, start_ops.elapsed().as_secs_f32());
    res
}

/// a^b mod p
#[time("debug", "Pow mod")]
pub fn pow_mod<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let mut res = server_key.create_trivial_radix(1, NB);
    let mut base = a.clone();
    let mut exponent = b.clone();
    let loop_end = <P as Numeric>::BITS;
    for _i in 0..loop_end {
        let _tmr = timer!(Level::Trace; "Pow Mod", "Bit {}", _i);

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
    }
    res
}

pub fn inverse_mod_pow<const NB: usize, P: Numeral>(
    a: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    let mut res = server_key.create_trivial_radix(1, NB);
    let mut base = a.clone();
    let mut exponent = p - P::cast_from(2u64);
    let loop_end = <P as Numeric>::BITS;
    for _i in 0..loop_end {
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
    }
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
        helper::set_client_key,
        numeral::Numeral,
        ops::{
            add_mod, double_mod, inverse_mod, inverse_mod_binary_gcd, inverse_mods,
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
        set_client_key(&client_key);
        const NUM_BLOCK: usize = 4;
        let p: u8 = 157;
        let a: u8 = 8;
        let b: u8 = 6;
        let e: u8 = 45;
        let f: u8 = 123;
        let h: u8 = 127;
        let i: u8 = 156;

        let c = inverse_mod_native(a, p);
        let enc_c =
            inverse_mod::<NUM_BLOCK, _>(&client_key.encrypt_radix(a, NUM_BLOCK), p, &server_key);
        assert_eq!(c, client_key.decrypt_radix::<u8>(&enc_c));

        let d = inverse_mod_native(b, p);
        let enc_d =
            inverse_mod::<NUM_BLOCK, _>(&client_key.encrypt_radix(b, NUM_BLOCK), p, &server_key);
        assert_eq!(d, client_key.decrypt_radix::<u8>(&enc_d));

        let m = inverse_mod_native(e, p);
        let enc_m =
            inverse_mod::<NUM_BLOCK, _>(&client_key.encrypt_radix(e, NUM_BLOCK), p, &server_key);
        assert_eq!(m, client_key.decrypt_radix::<u8>(&enc_m));

        let g = inverse_mod_native(f, p);
        let enc_g =
            inverse_mod::<NUM_BLOCK, _>(&client_key.encrypt_radix(f, NUM_BLOCK), p, &server_key);
        assert_eq!(g, client_key.decrypt_radix::<u8>(&enc_g));

        let j = inverse_mod_native(h, p);
        let enc_j =
            inverse_mod::<NUM_BLOCK, _>(&client_key.encrypt_radix(h, NUM_BLOCK), p, &server_key);
        assert_eq!(j, client_key.decrypt_radix::<u8>(&enc_j));

        let k = inverse_mod_native(i, p);
        let enc_k =
            inverse_mod::<NUM_BLOCK, _>(&client_key.encrypt_radix(i, NUM_BLOCK), p, &server_key);
        assert_eq!(k, client_key.decrypt_radix::<u8>(&enc_k));

        let l = inverse_mod_native(f, p);
        let enc_l =
            inverse_mod::<NUM_BLOCK, _>(&client_key.encrypt_radix(f, NUM_BLOCK), p, &server_key);
        assert_eq!(l, client_key.decrypt_radix::<u8>(&enc_l));
    }

    #[test]
    #[ignore = "bench"]
    fn bench_inverse() {
        let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);
        set_client_key(&client_key);
        const NUM_BLOCK: usize = 4;
        type Integer = u8;
        let p: Integer = 251;
        let a: Integer = 147;

        //const NUM_BLOCK: usize = 8;
        //type Integer = u16;
        //let p: Integer = 65521;
        //let a: Integer = 50725;

        let c = inverse_mod_native(a, p);

        let timer = Instant::now();
        let enc_c = inverse_mod_binary_gcd::<NUM_BLOCK, _>(
            &client_key.encrypt_radix(a, NUM_BLOCK),
            p,
            &server_key,
        );
        println!(
            "inverse_mod_binary_gcd: {:.2}s",
            timer.elapsed().as_secs_f32()
        );
        assert_eq!(c, client_key.decrypt_radix::<Integer>(&enc_c));

        let timer = Instant::now();
        let enc_c_other =
            inverse_mod::<NUM_BLOCK, _>(&client_key.encrypt_radix(a, NUM_BLOCK), p, &server_key);
        println!("inverse_mod: {:.2}s", timer.elapsed().as_secs_f32());
        assert_eq!(c, client_key.decrypt_radix::<Integer>(&enc_c_other));
    }

    #[test]
    fn correct_inverse_mods() {
        let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);
        set_client_key(&client_key);
        const NUM_BLOCK: usize = 4;
        let p: u8 = 157;
        let a: u8 = 8;
        let b: u8 = 123;

        let enc_a = client_key.encrypt_radix(a, NUM_BLOCK);
        let enc_b = client_key.encrypt_radix(b, NUM_BLOCK);

        let r_a = inverse_mod_native(a, p);
        let r_b = inverse_mod_native(b, p);

        let results = inverse_mods::<NUM_BLOCK, _>(&[enc_a, enc_b], p, &server_key);
        assert_eq!(r_a, client_key.decrypt_radix::<u8>(&results[0]));
        assert_eq!(r_b, client_key.decrypt_radix::<u8>(&results[1]));

        let results =
            inverse_mods::<NUM_BLOCK, _>(&[results[0].clone(), results[1].clone()], p, &server_key);
        assert_eq!(a, client_key.decrypt_radix::<u8>(&results[0]));
        assert_eq!(b, client_key.decrypt_radix::<u8>(&results[1]));
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
