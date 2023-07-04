use std::time::Instant;

use tfhe::integer::{
    block_decomposition::{DecomposableInto, RecomposableFrom},
    ClientKey, RadixCiphertext, ServerKey,
};

use crate::{
    helper::format,
    ops::{add_mod, double_mod, mul_mod, square_mod, sub_mod},
};

use super::inverse_mod;

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
    let (x_prime, z_prime) = rayon::join(
        || mul_mod::<NB, _>(&u, &w, p, server_key),
        || mul_mod::<NB, _>(&u, &square_mod::<NB, _>(&u, p, server_key), p, server_key),
    );
    let y_prime = sub_mod::<NB, _>(
        &mul_mod::<NB, _>(&t, &sub_mod::<NB, _>(&v, &w, p, server_key), p, server_key),
        &double_mod::<NB, _>(
            &square_mod::<NB, _>(&mul_mod::<NB, _>(&u, y, p, server_key), p, server_key),
            p,
            server_key,
        ),
        p,
        server_key,
    );

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
    let (t0, t1) = rayon::join(
        || mul_mod::<NB, _>(y0, z1, p, server_key),
        || mul_mod::<NB, _>(y1, z0, p, server_key),
    );
    // u0 = x0 * z1
    // u1 = x1 * z0
    let (u0, u1) = rayon::join(
        || mul_mod::<NB, _>(x0, z1, p, server_key),
        || mul_mod::<NB, _>(x1, z0, p, server_key),
    );
    // u = u0 - u1
    // t = t0 - t1
    let (u, t) = rayon::join(
        || sub_mod::<NB, _>(&u0, &u1, p, server_key),
        || sub_mod::<NB, _>(&t0, &t1, p, server_key),
    );
    // u2 = u^2
    // v = z0 * z1
    let (u2, v) = rayon::join(
        || square_mod::<NB, _>(&u, p, server_key),
        || mul_mod::<NB, _>(z0, z1, p, server_key),
    );
    // w = t^2 * v - u2 * (u0 + u1)
    let (t2v, u2u0u1) = rayon::join(
        || mul_mod::<NB, _>(&square_mod::<NB, _>(&t, p, server_key), &v, p, server_key),
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
    let ((u3, x_prime), tu2u0w) = rayon::join(
        || {
            rayon::join(
                || mul_mod::<NB, _>(&u, &u2, p, server_key),
                || mul_mod::<NB, _>(&u, &w, p, server_key),
            )
        },
        || {
            mul_mod::<NB, _>(
                &t,
                &sub_mod::<NB, _>(
                    &mul_mod::<NB, _>(&u2, &u0, p, server_key),
                    &w,
                    p,
                    server_key,
                ),
                p,
                server_key,
            )
        },
    );
    let (y_prime, z_prime) = rayon::join(
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

    (x_prime, y_prime, z_prime)
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
    let z_inv = inverse_mod::<NB, _>(z, p, server_key);

    rayon::join(
        || mul_mod::<NB, _>(&x, &z_inv, p, server_key),
        || mul_mod::<NB, _>(&y, &z_inv, p, server_key),
    )
}
