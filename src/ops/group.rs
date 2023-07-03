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
    P: DecomposableInto<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    x: &RadixCiphertext,
    y: &RadixCiphertext,
    z: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext, RadixCiphertext) {
    let timer = Instant::now();
    // y^2 = x^3+ax+b
    // in case a == 0

    // a = x^2
    // b = y^2
    let (a, b) = rayon::join(
        || square_mod::<NB, _>(x, p, server_key),
        || square_mod::<NB, _>(y, p, server_key),
    );
    println!("a,b - {:.2}s", timer.elapsed().as_secs_f64());
    // c = b^2
    let c = square_mod::<NB, _>(&b, p, server_key);
    println!("c - {:.2}s", timer.elapsed().as_secs_f64());
    // d = 2(2(x + b) - (a + c))
    let (xb2, ac) = rayon::join(
        || double_mod::<NB, _>(&add_mod::<NB, _>(x, &b, p, server_key), p, server_key),
        || add_mod::<NB, _>(&a, &c, p, server_key),
    );
    let d = double_mod::<NB, _>(
        &sub_mod::<NB, _>(
            &xb2, // 2(x + b)
            &ac,  // (a + c)
            p, server_key,
        ),
        p,
        server_key,
    );
    println!("d - {:.2}s", timer.elapsed().as_secs_f64());
    // e = 3a
    let e = add_mod::<NB, _>(&double_mod::<NB, _>(&a, p, server_key), &a, p, server_key);
    println!("e - {:.2}s", timer.elapsed().as_secs_f64());
    // f = e^2
    let f = square_mod::<NB, _>(&e, p, server_key);
    println!("f - {:.2}s", timer.elapsed().as_secs_f64());
    // z' = 2 * y1 * z1
    // x' = f - 2d
    let (z_prime, x_prime) = rayon::join(
        || double_mod::<NB, _>(&mul_mod::<NB, _>(y, z, p, server_key), p, server_key),
        || sub_mod::<NB, _>(&f, &double_mod::<NB, _>(&d, p, server_key), p, server_key),
    );
    println!("z',x' - {:.2}s", timer.elapsed().as_secs_f64());
    // y' = e(d - x') - 8c
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
            let enc_2c = double_mod::<NB, _>(&c, p, server_key);
            let enc_4c = double_mod::<NB, _>(&enc_2c, p, server_key);
            double_mod::<NB, _>(&enc_4c, p, server_key)
        },
    );
    let y_prime = sub_mod::<NB, _>(&edx, &c8, p, server_key);
    println!("y' - {:.2}s", timer.elapsed().as_secs_f64());

    (x_prime, y_prime, z_prime)
}

pub fn group_projective_add_affine<
    const NB: usize,
    P: DecomposableInto<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    x1: &RadixCiphertext,
    y1: &RadixCiphertext,
    z1: &RadixCiphertext,
    x2: &RadixCiphertext,
    y2: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext, RadixCiphertext) {
    // z1z1 = z1^2
    // y2z1 = y2 * z1
    let (z1z1, y2z1) = rayon::join(
        || square_mod::<NB, _>(z1, p, server_key),
        || mul_mod::<NB, _>(y2, z1, p, server_key),
    );
    // u2 = x2 * z1z1
    // s2 = y2z1 * z1z1
    let (u2, s2) = rayon::join(
        || mul_mod::<NB, _>(x2, &z1z1, p, server_key),
        || mul_mod::<NB, _>(&y2z1, &z1z1, p, server_key),
    );
    // h = u2 - x1
    let h = sub_mod::<NB, _>(&u2, x1, p, server_key);
    // hh = h^2
    let hh = square_mod::<NB, _>(&h, p, server_key);
    // i = 4 * hh
    let i = double_mod::<NB, _>(&double_mod::<NB, _>(&hh, p, server_key), p, server_key);
    // j = h * i
    let j = mul_mod::<NB, _>(&h, &i, p, server_key);
    // r = 2 * (s2 - y1)
    let r = double_mod::<NB, _>(&sub_mod::<NB, _>(&s2, y1, p, server_key), p, server_key);
    // v = x1 * i
    let v = mul_mod::<NB, _>(x1, &i, p, server_key);
    // x3 = r^2 - j - 2 * v
    let (r2, j2v) = rayon::join(
        || square_mod::<NB, _>(&r, p, server_key),
        || {
            let j2 = double_mod::<NB, _>(&j, p, server_key);
            sub_mod::<NB, _>(&j2, &double_mod::<NB, _>(&v, p, server_key), p, server_key)
        },
    );
    let x3 = sub_mod::<NB, _>(&r2, &j2v, p, server_key);
    // y3 = r * (v - x3) - 2 * y1 * j
    let (ryx3, y1j2) = rayon::join(
        || mul_mod::<NB, _>(&r, &sub_mod::<NB, _>(&v, &x3, p, server_key), p, server_key),
        || double_mod::<NB, _>(&mul_mod::<NB, _>(y1, &j, p, server_key), p, server_key),
    );
    let y3 = sub_mod::<NB, _>(&ryx3, &y1j2, p, server_key);
    // z3 = (z1 + h)^2 - (z1z1 + hh)
    let (z1h2, z1z1hh) = rayon::join(
        || square_mod::<NB, _>(&add_mod::<NB, _>(z1, &h, p, server_key), p, server_key),
        || add_mod::<NB, _>(z1, &h, p, server_key),
    );
    let z3 = sub_mod::<NB, _>(&z1h2, &z1z1hh, p, server_key);

    (x3, y3, z3)
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
    let z_inv2 = square_mod::<NB, _>(&z_inv, p, server_key);
    let z_inv3 = mul_mod::<NB, _>(&z_inv2, &z_inv, p, server_key);
    rayon::join(
        || mul_mod::<NB, _>(&x, &z_inv2, p, server_key),
        || mul_mod::<NB, _>(&y, &z_inv3, p, server_key),
    )
}
