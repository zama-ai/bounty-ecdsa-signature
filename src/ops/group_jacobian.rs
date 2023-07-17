use std::time::Instant;

use rand::Rng;
use rayon::prelude::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use tfhe::{
    core_crypto::prelude::Numeric,
    integer::{
        block_decomposition::{DecomposableInto, RecomposableFrom},
        RadixCiphertext, ServerKey,
    },
};

use crate::{
    helper::{format, read_client_key},
    ops::native::{add_mod_native, double_mod_native, mul_mod_native, sub_mod_native},
};

use super::{
    add_mod, double_mod, inverse_mod, mul_mod,
    native::{inverse_mod_native, square_mod_native},
    square_mod, sub_mod,
};

pub fn group_projective_double_native<
    P: DecomposableInto<u64> + RecomposableFrom<u8> + DecomposableInto<u8> + Copy + Sync,
>(
    x: P,
    y: P,
    z: P,
    p: P,
) -> (P, P, P) {
    // case curve a = 0
    // a = x^2
    let a = square_mod_native(x, p);
    // b = y^2
    let b = square_mod_native(y, p);
    // c = b^2
    let c = square_mod_native(b, p);
    // d = 2*((x + b)^2-(a + c))
    let xb2 = square_mod_native(add_mod_native(x, b, p), p);
    let ac = add_mod_native(a, c, p);
    let d = double_mod_native(sub_mod_native(xb2, ac, p), p);
    // e = 3*a
    let e = add_mod_native(double_mod_native(a, p), a, p);
    // f = e^2
    let f = square_mod_native(e, p);
    // z' = 2*y*z
    let z_prime = double_mod_native(mul_mod_native(y, z, p), p);
    // x' = f - 2*d
    let x_prime = sub_mod_native(f, double_mod_native(d, p), p);
    // y' = e*(d - x') - 8*c
    let edx = mul_mod_native(e, sub_mod_native(d, x_prime, p), p);
    let c2 = double_mod_native(c, p);
    let c4 = double_mod_native(c2, p);
    let c8 = double_mod_native(c4, p);
    let y_prime = sub_mod_native(edx, c8, p);

    (x_prime, y_prime, z_prime)
}

pub fn group_projective_add_affine_native<
    P: DecomposableInto<u64> + RecomposableFrom<u8> + DecomposableInto<u8> + Copy + Sync,
>(
    x: P,
    y: P,
    z: P,
    other_x: P,
    other_y: P,
    p: P,
) -> (P, P, P) {
    // z1z1 = z1^2
    let z1z1 = square_mod_native(z, p);
    // u2 = x2*z1z1
    let u2 = mul_mod_native(other_x, z1z1, p);
    // s2 = y2*z1*z1*z1
    let s2 = mul_mod_native(other_y, mul_mod_native(z1z1, z, p), p);

    if x == u2 && y == s2 {
        return group_projective_double_native(x, y, z, p);
    }

    // h = u2 - x1
    let h = sub_mod_native(u2, x, p);
    // hh = h^2
    let hh = square_mod_native(h, p);
    // i = 4*hh
    let i = double_mod_native(hh, p);
    let i = double_mod_native(i, p);
    // j = h*i
    let j = mul_mod_native(h, i, p);
    // r = 2*(s2 - y1)
    let r = double_mod_native(sub_mod_native(s2, y, p), p);
    // v = x1*i
    let v = mul_mod_native(x, i, p);
    // x3 = r^2 - j - 2*v
    let x3 = sub_mod_native(
        square_mod_native(r, p),
        add_mod_native(j, double_mod_native(v, p), p),
        p,
    );
    // y3 = r*(v - x3) - 2*y1*j
    let y3 = sub_mod_native(
        mul_mod_native(r, sub_mod_native(v, x3, p), p),
        double_mod_native(mul_mod_native(y, j, p), p),
        p,
    );
    // z3 = 2*z1*h
    let z3 = double_mod_native(mul_mod_native(z, h, p), p);

    (x3, y3, z3)
}

#[allow(clippy::too_many_arguments)]
pub fn group_projective_add_affine<
    const NB: usize,
    P: DecomposableInto<u64>
        + RecomposableFrom<u64>
        + RecomposableFrom<u8>
        + DecomposableInto<u8>
        + Copy
        + Sync,
>(
    x: &RadixCiphertext,
    y: &RadixCiphertext,
    z: &RadixCiphertext,
    other_x: &RadixCiphertext,
    other_y: &RadixCiphertext,
    other_flag_bit: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext, RadixCiphertext) {
    #[cfg(feature = "high_level_timing")]
    let ops_start = Instant::now();
    #[cfg(feature = "high_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);

    #[cfg(feature = "high_level_timing")]
    println!("group projective add affine start -- ref {}", task_ref);

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
    let ((mut x3, mut z3), yj2) = rayon::join(
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
    let mut y3 = sub_mod::<NB, _>(
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
        || server_key.smart_scalar_ne_parallelized(&mut z.clone(), 0),
        || server_key.smart_scalar_ne_parallelized(&mut other_flag_bit.clone(), 0),
    );
    server_key.trim_radix_blocks_msb_assign(&mut is_z0_non_zero, NB - 1);
    server_key.trim_radix_blocks_msb_assign(&mut is_z1_non_zero, NB - 1);
    let mut is_z0_z1_non_zero =
        server_key.smart_bitand_parallelized(&mut is_z0_non_zero, &mut is_z1_non_zero);
    let not_is_z0_z1_non_zero = server_key.smart_sub_parallelized(
        &mut server_key.create_trivial_radix(1, 1),
        &mut is_z0_z1_non_zero,
    );

    let (((mut xp1, mut xp2), (mut yp1, mut yp2)), (mut zp1, mut zp2)) = rayon::join(
        || {
            rayon::join(
                || {
                    rayon::join(
                        || {
                            server_key
                                .smart_mul_parallelized(&mut x3, &mut is_z0_z1_non_zero.clone())
                        },
                        || {
                            server_key.smart_mul_parallelized(
                                &mut server_key
                                    .smart_add_parallelized(&mut x.clone(), &mut other_x.clone()),
                                &mut not_is_z0_z1_non_zero.clone(),
                            )
                        },
                    )
                },
                || {
                    rayon::join(
                        || {
                            server_key
                                .smart_mul_parallelized(&mut y3, &mut is_z0_z1_non_zero.clone())
                        },
                        || {
                            server_key.smart_mul_parallelized(
                                &mut server_key
                                    .smart_add_parallelized(&mut y.clone(), &mut other_y.clone()),
                                &mut not_is_z0_z1_non_zero.clone(),
                            )
                        },
                    )
                },
            )
        },
        || {
            rayon::join(
                || server_key.smart_mul_parallelized(&mut z3, &mut is_z0_z1_non_zero.clone()),
                || {
                    server_key.smart_mul_parallelized(
                        &mut server_key
                            .smart_add_parallelized(&mut z.clone(), &mut other_flag_bit.clone()),
                        &mut not_is_z0_z1_non_zero.clone(),
                    )
                },
            )
        },
    );

    let ((x_prime, y_prime), z_prime) = rayon::join(
        || {
            rayon::join(
                || server_key.smart_add_parallelized(&mut xp1, &mut xp2),
                || server_key.smart_add_parallelized(&mut yp1, &mut yp2),
            )
        },
        || server_key.smart_add_parallelized(&mut zp1, &mut zp2),
    );

    #[cfg(feature = "high_level_timing")]
    println!(
        "group projective add affine done in {:.2}s -- ref {}",
        ops_start.elapsed().as_secs_f64(),
        task_ref
    );

    (x_prime, y_prime, z_prime)
}

pub fn group_projective_double<
    const NB: usize,
    P: DecomposableInto<u64> + RecomposableFrom<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    x: &RadixCiphertext,
    y: &RadixCiphertext,
    z: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext, RadixCiphertext) {
    // case curve a = 0
    // a = x^2
    // b = y^2

    #[cfg(feature = "high_level_timing")]
    let ops_start = Instant::now();
    #[cfg(feature = "high_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);

    #[cfg(feature = "high_level_timing")]
    println!("group projective double jacobian start -- ref {}", task_ref);

    let (a, b) = rayon::join(
        || square_mod::<NB, _>(x, p, server_key),
        || square_mod::<NB, _>(y, p, server_key),
    );
    // c = b^2
    let c = square_mod::<NB, _>(&b, p, server_key);
    // d = 2*((x + b)^2-(a + c))
    let (xb, ac) = rayon::join(
        || add_mod::<NB, _>(x, &b, p, server_key),
        || add_mod::<NB, _>(&a, &c, p, server_key),
    );
    let d = double_mod::<NB, _>(
        &sub_mod::<NB, _>(&square_mod::<NB, _>(&xb, p, server_key), &ac, p, server_key),
        p,
        server_key,
    );
    // e = 3*a
    let e = add_mod::<NB, _>(&double_mod::<NB, _>(&a, p, server_key), &a, p, server_key);
    // f = e^2
    // z' = 2*y*z
    let (f, z_prime) = rayon::join(
        || square_mod::<NB, _>(&e, p, server_key),
        || double_mod::<NB, _>(&mul_mod::<NB, _>(y, z, p, server_key), p, server_key),
    );
    // x' = f - 2*d
    let x_prime = sub_mod::<NB, _>(&f, &double_mod::<NB, _>(&d, p, server_key), p, server_key);
    // y' = e*(d - x') - 8*c
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
            let c2 = double_mod::<NB, _>(&c, p, server_key);
            let c4 = double_mod::<NB, _>(&c2, p, server_key);
            double_mod::<NB, _>(&c4, p, server_key)
        },
    );
    let y_prime = sub_mod::<NB, _>(&edx, &c8, p, server_key);
    #[cfg(feature = "high_level_timing")]
    println!(
        "group projective double done in {:.2}s -- ref {}",
        ops_start.elapsed().as_secs_f64(),
        task_ref
    );
    (x_prime, y_prime, z_prime)
}

#[allow(clippy::too_many_arguments)]
pub fn group_projective_add_projective<
    const NB: usize,
    P: DecomposableInto<u64> + RecomposableFrom<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    x0: &RadixCiphertext,
    y0: &RadixCiphertext,
    z0: &RadixCiphertext,
    x1: &RadixCiphertext,
    y1: &RadixCiphertext,
    z1: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext, RadixCiphertext) {
    #[cfg(feature = "high_level_timing")]
    let ops_start = Instant::now();
    #[cfg(feature = "high_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);
    #[cfg(feature = "high_level_timing")]
    println!("group projective add jacobian start -- ref {}", task_ref);

    // z0z0 = z0^2
    // z1z1 = z1^2
    let (z0z0, z1z1) = rayon::join(
        || square_mod::<NB, _>(z0, p, server_key),
        || square_mod::<NB, _>(z1, p, server_key),
    );
    // u0 = x0*z1z1
    // u1 = x1*z0z0
    let (u0, u1) = rayon::join(
        || mul_mod::<NB, _>(x0, &z1z1, p, server_key),
        || mul_mod::<NB, _>(x1, &z0z0, p, server_key),
    );
    // s0 = y0*z1*z1z1
    // s1 = y1*z0*z0z0
    let (s0, s1) = rayon::join(
        || {
            mul_mod::<NB, _>(
                y0,
                &mul_mod::<NB, _>(z1, &z1z1, p, server_key),
                p,
                server_key,
            )
        },
        || {
            mul_mod::<NB, _>(
                y1,
                &mul_mod::<NB, _>(z0, &z0z0, p, server_key),
                p,
                server_key,
            )
        },
    );
    // h = u1 - u0
    // r = 2*(s1 - s0)
    let (h, r) = rayon::join(
        || sub_mod::<NB, _>(&u1, &u0, p, server_key),
        || double_mod::<NB, _>(&sub_mod::<NB, _>(&s1, &s0, p, server_key), p, server_key),
    );
    // i = (2*h)^2
    let i = square_mod::<NB, _>(&double_mod::<NB, _>(&h, p, server_key), p, server_key);
    // j = h*i
    // v = u0*i
    let (j, v) = rayon::join(
        || mul_mod::<NB, _>(&h, &i, p, server_key),
        || mul_mod::<NB, _>(&u0, &i, p, server_key),
    );
    // x_prime = r^2 - j - 2*v
    // y_prime = r*(v - x_prime) - 2*s0*j
    // z_prime = ((z0 + z1)^2 - z0z0 - z1z1)*h
    let ((r2, s0j2), z0z12) = rayon::join(
        || {
            rayon::join(
                || square_mod::<NB, _>(&r, p, server_key),
                || double_mod::<NB, _>(&mul_mod::<NB, _>(&s0, &j, p, server_key), p, server_key),
            )
        },
        || square_mod::<NB, _>(&add_mod::<NB, _>(z0, z1, p, server_key), p, server_key),
    );
    let mut x_prime = sub_mod::<NB, _>(
        &sub_mod::<NB, _>(&r2, &j, p, server_key),
        &double_mod::<NB, _>(&v, p, server_key),
        p,
        server_key,
    );
    let (mut y_prime, mut z_prime) = rayon::join(
        || {
            sub_mod::<NB, _>(
                &mul_mod::<NB, _>(
                    &r,
                    &sub_mod::<NB, _>(&v, &x_prime, p, server_key),
                    p,
                    server_key,
                ),
                &s0j2,
                p,
                server_key,
            )
        },
        || {
            mul_mod::<NB, _>(
                &sub_mod::<NB, _>(
                    &sub_mod::<NB, _>(&z0z12, &z0z0, p, server_key),
                    &z1z1,
                    p,
                    server_key,
                ),
                &h,
                p,
                server_key,
            )
        },
    );

    // z1'/z0' 0  1
    //    0    x' x1
    //    1    x0 x0
    // x'' =  x' * is_z0_z1_non_zero + (x0 + x1) * not_is_z0_z1_non_zero
    // y'' =  y' * is_z0_z1_non_zero + (y0 + y1) * not_is_z0_z1_non_zero
    // z'' =  z' * is_z0_z1_non_zero + (z0 + z1) * not_is_z0_z1_non_zero
    let (mut is_z0_non_zero, mut is_z1_non_zero) = rayon::join(
        || server_key.smart_scalar_ne_parallelized(&mut z0.clone(), 0),
        || server_key.smart_scalar_ne_parallelized(&mut z1.clone(), 0),
    );
    server_key.trim_radix_blocks_msb_assign(&mut is_z0_non_zero, NB - 1);
    server_key.trim_radix_blocks_msb_assign(&mut is_z1_non_zero, NB - 1);
    let mut is_z0_z1_non_zero =
        server_key.smart_bitand_parallelized(&mut is_z0_non_zero, &mut is_z1_non_zero);
    let not_is_z0_z1_non_zero = server_key.smart_sub_parallelized(
        &mut server_key.create_trivial_radix(1, 1),
        &mut is_z0_z1_non_zero,
    );

    let (((mut xp1, mut xp2), (mut yp1, mut yp2)), (mut zp1, mut zp2)) = rayon::join(
        || {
            rayon::join(
                || {
                    rayon::join(
                        || {
                            server_key.smart_mul_parallelized(
                                &mut x_prime,
                                &mut is_z0_z1_non_zero.clone(),
                            )
                        },
                        || {
                            server_key.smart_mul_parallelized(
                                &mut server_key
                                    .smart_add_parallelized(&mut x0.clone(), &mut x1.clone()),
                                &mut not_is_z0_z1_non_zero.clone(),
                            )
                        },
                    )
                },
                || {
                    rayon::join(
                        || {
                            server_key.smart_mul_parallelized(
                                &mut y_prime,
                                &mut is_z0_z1_non_zero.clone(),
                            )
                        },
                        || {
                            server_key.smart_mul_parallelized(
                                &mut server_key
                                    .smart_add_parallelized(&mut y0.clone(), &mut y1.clone()),
                                &mut not_is_z0_z1_non_zero.clone(),
                            )
                        },
                    )
                },
            )
        },
        || {
            rayon::join(
                || server_key.smart_mul_parallelized(&mut z_prime, &mut is_z0_z1_non_zero.clone()),
                || {
                    server_key.smart_mul_parallelized(
                        &mut server_key.smart_add_parallelized(&mut z0.clone(), &mut z1.clone()),
                        &mut not_is_z0_z1_non_zero.clone(),
                    )
                },
            )
        },
    );

    ((x_prime, y_prime), z_prime) = rayon::join(
        || {
            rayon::join(
                || server_key.smart_add_parallelized(&mut xp1, &mut xp2),
                || server_key.smart_add_parallelized(&mut yp1, &mut yp2),
            )
        },
        || server_key.smart_add_parallelized(&mut zp1, &mut zp2),
    );

    #[cfg(feature = "high_level_timing")]
    println!(
        "group projective add done in {:.2}s -- ref {}",
        ops_start.elapsed().as_secs_f64(),
        task_ref
    );

    (x_prime, y_prime, z_prime)
}

pub fn group_projective_scalar_mul<
    const NB: usize,
    P: DecomposableInto<u64> + RecomposableFrom<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    x: &RadixCiphertext,
    y: &RadixCiphertext,
    z: &RadixCiphertext,
    scalar: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext, RadixCiphertext) {
    #[cfg(feature = "high_level_timing")]
    let ops_start = Instant::now();
    #[cfg(feature = "high_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);
    #[cfg(feature = "high_level_timing")]
    println!(
        "group projective scalar mul jacobian start -- ref {}",
        task_ref
    );

    let mut tmp_x = x.clone();
    let mut tmp_y = y.clone();
    let mut tmp_z = z.clone();
    let mut scalar = scalar.clone();
    let mut res_x = server_key.create_trivial_radix(0, NB);
    let mut res_y = server_key.create_trivial_radix(0, NB);
    let mut res_z = server_key.create_trivial_radix(0, NB);

    for _i in 0..<P as Numeric>::BITS {
        #[cfg(feature = "high_level_timing")]
        let bit_start = Instant::now();

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
                            || {
                                server_key
                                    .smart_mul_parallelized(&mut tmp_x.clone(), &mut bit.clone())
                            },
                            || {
                                server_key
                                    .smart_mul_parallelized(&mut tmp_y.clone(), &mut bit.clone())
                            },
                        )
                    },
                    || server_key.smart_mul_parallelized(&mut tmp_z.clone(), &mut bit.clone()),
                );
                group_projective_add_projective::<NB, _>(
                    &res_x, &res_y, &res_z, &x_to_add, &y_to_add, &z_to_add, p, server_key,
                )
            },
            || group_projective_double::<NB, _>(&tmp_x, &tmp_y, &tmp_z, p, server_key),
        );
        #[cfg(feature = "high_level_timing")]
        read_client_key(|client_key| {
            println!("Bit = {}", format(client_key.decrypt_radix::<P>(&bit)),);
            println!(
                "Res {},{},{}",
                format(client_key.decrypt_radix::<P>(&res_x)),
                format(client_key.decrypt_radix::<P>(&res_y)),
                format(client_key.decrypt_radix::<P>(&res_z)),
            );
            println!(
                "Tmp {},{},{}",
                format(client_key.decrypt_radix::<P>(&tmp_x)),
                format(client_key.decrypt_radix::<P>(&tmp_y)),
                format(client_key.decrypt_radix::<P>(&tmp_z)),
            );
            println!(
                "----Scalar mul bit {_i} done in {:.2}s -- ref {}",
                bit_start.elapsed().as_secs_f32(),
                task_ref
            );
        });
    }

    #[cfg(feature = "high_level_timing")]
    println!(
        "group projective scalar mul done in {:.2}s -- ref {}",
        ops_start.elapsed().as_secs_f64(),
        task_ref
    );
    (res_x, res_y, res_z)
}

pub fn group_projective_scalar_mul_constant<
    const NB: usize,
    P: DecomposableInto<u64>
        + RecomposableFrom<u64>
        + RecomposableFrom<u8>
        + DecomposableInto<u8>
        + Copy
        + Sync,
>(
    x: P,
    y: P,
    scalar: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext, RadixCiphertext) {
    #[cfg(feature = "high_level_timing")]
    let ops_start = Instant::now();
    #[cfg(feature = "high_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);
    #[cfg(feature = "high_level_timing")]
    println!(
        "group projective scalar mul jacobian start -- ref {}",
        task_ref
    );

    let mut tmp_x = x;
    let mut tmp_y = y;
    let mut scalar = scalar.clone();
    let mut res_x = server_key.create_trivial_radix(0, NB);
    let mut res_y = server_key.create_trivial_radix(0, NB);
    let mut res_z = server_key.create_trivial_radix(0, NB);

    for _i in 0..<P as Numeric>::BITS {
        #[cfg(feature = "high_level_timing")]
        let bit_start = Instant::now();

        let (mut bit, new_scalar) = rayon::join(
            || server_key.scalar_bitand_parallelized(&scalar, 1),
            || server_key.scalar_right_shift_parallelized(&scalar, 1),
        );
        server_key.trim_radix_blocks_msb_assign(&mut bit, NB - 1);
        scalar = new_scalar;

        let (x_to_add, y_to_add) = rayon::join(
            || {
                server_key.smart_mul_parallelized(
                    &mut server_key.create_trivial_radix(tmp_x, NB),
                    &mut bit.clone(),
                )
            },
            || {
                server_key.smart_mul_parallelized(
                    &mut server_key.create_trivial_radix(tmp_y, NB),
                    &mut bit.clone(),
                )
            },
        );

        (res_x, res_y, res_z) = group_projective_add_affine::<NB, _>(
            &res_x, &res_y, &res_z, &x_to_add, &y_to_add, &bit, p, server_key,
        );

        (tmp_x, tmp_y) = {
            let (tmp_x_new, temp_y_new, temp_z_new) =
                group_projective_double_native(tmp_x, tmp_y, P::ONE, p);
            group_projective_into_affine_native(tmp_x_new, temp_y_new, temp_z_new, p)
        };

        #[cfg(feature = "high_level_timing")]
        read_client_key(|client_key| {
            println!("Bit = {}", format(client_key.decrypt_radix::<P>(&bit)),);
            println!(
                "Res {},{},{}",
                format(client_key.decrypt_radix::<P>(&res_x)),
                format(client_key.decrypt_radix::<P>(&res_y)),
                format(client_key.decrypt_radix::<P>(&res_z)),
            );
            println!("Tmp {},{}", format(tmp_x), format(tmp_y),);
            println!(
                "----Scalar mul bit {_i} done in {:.2}s -- ref {}",
                bit_start.elapsed().as_secs_f32(),
                task_ref
            );
        });
    }

    #[cfg(feature = "high_level_timing")]
    println!(
        "group projective scalar mul done in {:.2}s -- ref {}",
        ops_start.elapsed().as_secs_f64(),
        task_ref
    );
    (res_x, res_y, res_z)
}

pub fn group_projective_scalar_mul_constant_windowed<
    const W: usize,
    const NB: usize,
    P: DecomposableInto<u64>
        + RecomposableFrom<u64>
        + RecomposableFrom<u8>
        + DecomposableInto<u8>
        + Copy
        + Sync,
>(
    x: P,
    y: P,
    scalar: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext, RadixCiphertext) {
    #[cfg(feature = "high_level_timing")]
    let ops_start = Instant::now();
    #[cfg(feature = "high_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);
    #[cfg(feature = "high_level_timing")]
    println!(
        "group projective scalar mul jacobian start -- ref {}",
        task_ref
    );

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

        #[cfg(feature = "high_level_timing")]
        println!("----Start: Scalar mul bit {_ic:?}");
        #[cfg(feature = "high_level_timing")]
        let bit_start = Instant::now();

        // get the next W bits
        let mut tmp_bits = vec![
            (
                server_key.create_trivial_radix(0, NB),
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
                    server_key.smart_sub_parallelized(
                        &mut server_key.create_trivial_radix(P::ONE, 1),
                        &mut bit,
                    ),
                    bit,
                    shifted,
                )
            })
            .collect_into_vec(&mut tmp_bits);
        let mut bits = vec![];
        let mut not_bits = vec![];
        for (bit, not_bit, shifted) in tmp_bits {
            bits.push(bit);
            not_bits.push(not_bit);
            scalar = shifted;
        }
        #[cfg(feature = "high_level_timing")]
        println!(
            "----Calculating bits done in {:.2}s",
            bit_start.elapsed().as_secs_f32()
        );

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
        let mut selected_points = (
            server_key.create_trivial_radix(0, NB),
            server_key.create_trivial_radix(0, NB),
        );
        #[cfg(feature = "high_level_timing")]
        let now = Instant::now();
        for (i, point) in points
            .iter()
            .enumerate()
            .take(2usize.pow(chunk_size as u32))
            .skip(1)
        {
            #[cfg(feature = "high_level_timing")]
            let now = Instant::now();
            let mut selected_bit = match i & 1 == 0 {
                true => not_bits[0].clone(),
                false => bits[0].clone(),
            };
            for j in 1..chunk_size {
                let mut selected_bit_and = match i & 2usize.pow(j as u32) == 0 {
                    true => not_bits[j].clone(),
                    false => bits[j].clone(),
                };
                server_key
                    .smart_bitand_assign_parallelized(&mut selected_bit, &mut selected_bit_and);
            }
            rayon::join(
                || {
                    server_key.smart_add_assign_parallelized(
                        &mut selected_points.0,
                        &mut server_key.smart_mul_parallelized(
                            &mut server_key.create_trivial_radix(point.0, NB),
                            &mut selected_bit.clone(),
                        ),
                    )
                },
                || {
                    server_key.smart_add_assign_parallelized(
                        &mut selected_points.1,
                        &mut server_key.smart_mul_parallelized(
                            &mut server_key.create_trivial_radix(point.1, NB),
                            &mut selected_bit.clone(),
                        ),
                    )
                },
            );
            #[cfg(feature = "high_level_timing")]
            if i % 100 == 1 {
                println!(
                    "----Calculating selector option {i} done in {:.2}s",
                    now.elapsed().as_secs_f32()
                );
            }
        }
        #[cfg(feature = "high_level_timing")]
        println!(
            "----Calculating selector done in {:.2}s",
            now.elapsed().as_secs_f32()
        );

        // check if all bits are not zero for flag bit
        let mut all_not_zero = bits[0].clone();
        for bit in bits.iter_mut().take(chunk_size).skip(1) {
            server_key.smart_bitor_assign_parallelized(&mut all_not_zero, bit);
        }

        let (_res_x_before, _res_y_before, _res_z_before) =
            (res_x.clone(), res_y.clone(), res_z.clone());

        // add the point
        (res_x, res_y, res_z) = group_projective_add_affine::<NB, _>(
            &res_x,
            &res_y,
            &res_z,
            &selected_points.0,
            &selected_points.1,
            &all_not_zero,
            p,
            server_key,
        );

        #[cfg(feature = "high_level_timing")]
        read_client_key(|client_key| {
            println!(
                "Bits = {:?}",
                bits.iter()
                    .map(|bit| format(client_key.decrypt_radix::<P>(bit)))
                    .collect::<Vec<_>>()
            );

            println!(
                "Res {},{},{}",
                format(client_key.decrypt_radix::<P>(&res_x)),
                format(client_key.decrypt_radix::<P>(&res_y)),
                format(client_key.decrypt_radix::<P>(&res_z)),
            );
            println!(
                "Selected = {},{}",
                format(client_key.decrypt_radix::<P>(&selected_points.0)),
                format(client_key.decrypt_radix::<P>(&selected_points.1))
            );
            println!(
                "----Scalar mul bit {_ic:?} done in {:.2}s -- ref {}",
                bit_start.elapsed().as_secs_f32(),
                task_ref
            );
        });
    }

    #[cfg(feature = "high_level_timing")]
    println!(
        "group projective scalar mul done in {:.2}s -- ref {}",
        ops_start.elapsed().as_secs_f64(),
        task_ref
    );
    (res_x, res_y, res_z)
}

pub fn group_projective_into_affine<
    const NB: usize,
    P: DecomposableInto<u64> + RecomposableFrom<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    x: &RadixCiphertext,
    y: &RadixCiphertext,
    z: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext) {
    #[cfg(feature = "high_level_timing")]
    let ops_start = Instant::now();
    #[cfg(feature = "high_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);
    #[cfg(feature = "high_level_timing")]
    println!(
        "group projective into affine jacobian start -- ref {}",
        task_ref
    );

    let z_inv = inverse_mod::<NB, _>(z, p, server_key);
    let z_inv2 = square_mod::<NB, _>(&z_inv, p, server_key);
    let z_inv3 = mul_mod::<NB, _>(&z_inv2, &z_inv, p, server_key);

    let res = rayon::join(
        || mul_mod::<NB, _>(x, &z_inv2, p, server_key),
        || mul_mod::<NB, _>(y, &z_inv3, p, server_key),
    );
    #[cfg(feature = "high_level_timing")]
    println!(
        "group projective into affine done in {:.2}s -- ref {}",
        ops_start.elapsed().as_secs_f64(),
        task_ref
    );
    res
}

pub fn group_projective_into_affine_native<
    P: DecomposableInto<u64> + RecomposableFrom<u8> + DecomposableInto<u8> + Copy + Sync,
>(
    x: P,
    y: P,
    z: P,
    p: P,
) -> (P, P) {
    let z_inv = inverse_mod_native(z, p);
    let z_inv2 = square_mod_native(z_inv, p);
    let z_inv3 = mul_mod_native(z_inv2, z_inv, p);

    (mul_mod_native(x, z_inv2, p), mul_mod_native(y, z_inv3, p))
}

#[cfg(test)]
mod tests {

    use tfhe::{integer::keycache::IntegerKeyCache, shortint::prelude::PARAM_MESSAGE_2_CARRY_2};

    use crate::ops::group_jacobian::{
        group_projective_add_affine, group_projective_add_affine_native, group_projective_double,
        group_projective_double_native, group_projective_into_affine_native,
    };

    #[test]
    fn correct_jacobian_double() {
        let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);

        const NUM_BLOCK: usize = 4;
        type Integer = u8;
        let p: Integer = 251;
        let x1: Integer = 8;
        let y1: Integer = 45;

        let ct_x1 = client_key.encrypt_radix(x1, NUM_BLOCK);
        let ct_y1 = client_key.encrypt_radix(y1, NUM_BLOCK);

        let (x_new, y_new, z_new) = group_projective_double::<NUM_BLOCK, _>(
            &ct_x1,
            &ct_y1,
            &server_key.create_trivial_radix(1, NUM_BLOCK),
            p,
            &server_key,
        );
        let x_dec = client_key.decrypt_radix::<Integer>(&x_new);
        let y_dec = client_key.decrypt_radix::<Integer>(&y_new);
        let z_dec = client_key.decrypt_radix::<Integer>(&z_new);

        assert_eq!(x_dec, 134);
        assert_eq!(y_dec, 104);
        assert_eq!(z_dec, 90);
    }

    #[test]
    fn correct_jacobian_add_affine() {
        let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);

        const NUM_BLOCK: usize = 4;
        type Integer = u8;
        let p: Integer = 251;
        let x1: Integer = 48;
        let y1: Integer = 68;
        let z1: Integer = 153;
        let x2: Integer = 56;
        let y2: Integer = 225;

        let ct_x1 = client_key.encrypt_radix(x1, NUM_BLOCK);
        let ct_y1 = client_key.encrypt_radix(y1, NUM_BLOCK);
        let ct_z1 = client_key.encrypt_radix(z1, NUM_BLOCK);
        let ct_x2 = client_key.encrypt_radix(x2, NUM_BLOCK);
        let ct_y2 = client_key.encrypt_radix(y2, NUM_BLOCK);

        let (x_new, y_new, z_new) = group_projective_add_affine::<NUM_BLOCK, _>(
            &ct_x1,
            &ct_y1,
            &ct_z1,
            &ct_x2,
            &ct_y2,
            &client_key.encrypt_radix(1, NUM_BLOCK),
            p,
            &server_key,
        );
        let x_dec = client_key.decrypt_radix::<Integer>(&x_new);
        let y_dec = client_key.decrypt_radix::<Integer>(&y_new);
        let z_dec = client_key.decrypt_radix::<Integer>(&z_new);

        let res = group_projective_add_affine_native(x1, y1, z1, x2, y2, p);

        assert_eq!(x_dec, res.0);
        assert_eq!(y_dec, res.1);
        assert_eq!(z_dec, res.2);
    }

    #[test]
    fn correct_native_group_ops_jacobian() {
        let p: u8 = 251;
        let x: u8 = 8;
        let y: u8 = 45;
        let z: u8 = 1;

        let (xp, yp, zp) = group_projective_double_native(x, y, z, p);
        let (xn, yn) = group_projective_into_affine_native(xp, yp, zp, p);

        assert_eq!(xn, 157);
        assert_eq!(yn, 22);
    }
}
