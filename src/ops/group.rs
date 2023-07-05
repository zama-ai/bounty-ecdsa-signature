use std::time::Instant;

use tfhe::{
    core_crypto::prelude::Numeric,
    integer::{
        block_decomposition::{DecomposableInto, RecomposableFrom},
        ClientKey, RadixCiphertext, ServerKey,
    },
};

use crate::{
    helper::format,
    ops::{add_mod, double_mod, mul_mod, square_mod, sub_mod},
};

use super::inverse_mod;

/// Projective point at infinity in homogeneous coordinates
pub fn group_projective_zero_homogenous<const NB: usize>(
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext, RadixCiphertext) {
    (
        server_key.create_trivial_radix(0, NB),
        server_key.create_trivial_radix(1, NB),
        server_key.create_trivial_radix(0, NB),
    )
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
    // u = 2yz
    // t = 3x^2 + a * z^2 -> a = 0 so t = 3x^2
    let (u, t) = rayon::join(
        || double_mod::<NB, _>(&mul_mod::<NB, _>(y, z, p, server_key), p, server_key),
        || {
            let t = square_mod::<NB, _>(x, p, server_key);
            add_mod::<NB, _>(&double_mod::<NB, _>(&t, p, server_key), &t, p, server_key)
        },
    );
    // v = 2uxy
    let v = double_mod::<NB, _>(
        &mul_mod::<NB, _>(&u, &mul_mod::<NB, _>(x, y, p, server_key), p, server_key),
        p,
        server_key,
    );
    // w = t^2 - 2v
    let w = {
        let (t2, v2) = rayon::join(
            || square_mod::<NB, _>(&t, p, server_key),
            || double_mod::<NB, _>(&v, p, server_key),
        );
        sub_mod::<NB, _>(&t2, &v2, p, server_key)
    };
    // x' = uw
    // y' = t(v - w) - 2(uy)^2
    // z' = u^3
    let ((x_prime, u2), (tvw, uy)) = rayon::join(
        || {
            rayon::join(
                || mul_mod::<NB, _>(&u, &w, p, server_key),
                || square_mod::<NB, _>(&u, p, server_key),
            )
        },
        || {
            rayon::join(
                || mul_mod::<NB, _>(&t, &sub_mod::<NB, _>(&v, &w, p, server_key), p, server_key),
                || mul_mod::<NB, _>(&u, y, p, server_key),
            )
        },
    );
    let (uy2, z_prime) = rayon::join(
        || double_mod::<NB, _>(&square_mod::<NB, _>(&uy, p, server_key), p, server_key),
        || mul_mod::<NB, _>(&u2, &u, p, server_key),
    );
    let y_prime = sub_mod::<NB, _>(&tvw, &uy2, p, server_key);

    (x_prime, y_prime, z_prime)
}

pub fn group_projective_add_affine<
    const NB: usize,
    P: DecomposableInto<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    _x0: &RadixCiphertext,
    _y0: &RadixCiphertext,
    _z0: &RadixCiphertext,
    _x1: &RadixCiphertext,
    _y1: &RadixCiphertext,
    _z1: &RadixCiphertext,
    _p: P,
    _server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext, RadixCiphertext) {
    todo!()
}

pub fn group_projective_add_projective<
    const NB: usize,
    P: DecomposableInto<u64> + DecomposableInto<u8> + Copy + Sync,
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
    // t0 = y0 * z1
    // t1 = y1 * z0
    // u0 = x0 * z1
    // u1 = x1 * z0
    let ((t0, t1), (u0, u1)) = rayon::join(
        || {
            rayon::join(
                || mul_mod::<NB, _>(y0, z1, p, server_key),
                || mul_mod::<NB, _>(y1, z0, p, server_key),
            )
        },
        || {
            rayon::join(
                || mul_mod::<NB, _>(x0, z1, p, server_key),
                || mul_mod::<NB, _>(x1, z0, p, server_key),
            )
        },
    );
    // u = u0 - u1
    // t = t0 - t1
    let (u, t) = rayon::join(
        || sub_mod::<NB, _>(&u0, &u1, p, server_key),
        || sub_mod::<NB, _>(&t0, &t1, p, server_key),
    );
    // u2 = u^2
    // v = z0 * z1
    // w = t^2 * v - u2 * (u0 + u1)
    let ((u2, v), t2) = rayon::join(
        || {
            rayon::join(
                || square_mod::<NB, _>(&u, p, server_key),
                || mul_mod::<NB, _>(z0, z1, p, server_key),
            )
        },
        || square_mod::<NB, _>(&t, p, server_key),
    );
    let (t2v, u2u0u1) = rayon::join(
        || mul_mod::<NB, _>(&t2, &v, p, server_key),
        || {
            mul_mod::<NB, _>(
                &u2,
                &add_mod::<NB, _>(&u0, &u1, p, server_key),
                p,
                server_key,
            )
        },
    );
    let w = sub_mod::<NB, _>(&t2v, &u2u0u1, p, server_key);
    // u3 = u * u2
    // x' = u * w
    // y' = t * (u2 * u0 - w) - u3 * t0
    // z' = u3 * v
    let ((u3, mut x_prime), u2u0w) = rayon::join(
        || {
            rayon::join(
                || mul_mod::<NB, _>(&u, &u2, p, server_key),
                || mul_mod::<NB, _>(&u, &w, p, server_key),
            )
        },
        || {
            sub_mod::<NB, _>(
                &mul_mod::<NB, _>(&u2, &u0, p, server_key),
                &w,
                p,
                server_key,
            )
        },
    );
    let tu2u0w = mul_mod::<NB, _>(&t, &u2u0w, p, server_key);
    let (mut y_prime, mut z_prime) = rayon::join(
        || {
            sub_mod::<NB, _>(
                &tu2u0w,
                &mul_mod::<NB, _>(&u3, &t0, p, server_key),
                p,
                server_key,
            )
        },
        || mul_mod::<NB, _>(&u3, &v, p, server_key),
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

    (x_prime, y_prime, z_prime)
}

pub fn group_projective_into_affine<
    const NB: usize,
    P: DecomposableInto<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    x: &RadixCiphertext,
    y: &RadixCiphertext,
    z: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext) {
    let z_inv = inverse_mod::<NB, _>(z, p, server_key);

    rayon::join(
        || mul_mod::<NB, _>(&x, &z_inv, p, server_key),
        || mul_mod::<NB, _>(&y, &z_inv, p, server_key),
    )
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
    client_key: &ClientKey,
) -> (RadixCiphertext, RadixCiphertext, RadixCiphertext) {
    let mut tmp_x = x.clone();
    let mut tmp_y = y.clone();
    let mut tmp_z = z.clone();
    let mut scalar = scalar.clone();
    let mut res_x = server_key.create_trivial_radix(0, NB);
    let mut res_y = server_key.create_trivial_radix(0, NB);
    let mut res_z = server_key.create_trivial_radix(0, NB);

    for i in 0..<P as Numeric>::BITS {
        let now = Instant::now();
        let (mut bit, new_scalar) = rayon::join(
            || server_key.scalar_bitand_parallelized(&mut scalar.clone(), 1),
            || server_key.scalar_right_shift_parallelized(&mut scalar.clone(), 1),
        );
        server_key.trim_radix_blocks_msb(&mut bit, NB - 1);
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
            "Scalar mul bit {} took {:.2}s",
            i,
            now.elapsed().as_secs_f32()
        );
    }

    (res_x, res_y, res_z)
}
