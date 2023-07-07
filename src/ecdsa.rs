use std::time::Instant;

use tfhe::{
    core_crypto::prelude::Numeric,
    integer::{
        block_decomposition::{DecomposableInto, RecomposableFrom},
        ClientKey, RadixCiphertext, ServerKey,
    },
};

use crate::{
    helper::{format, read_client_key},
    ops::{
        add_mod,
        group::{
            group_projective_add_projective, group_projective_into_affine,
            group_projective_scalar_mul,
        },
        inverse_mod, modulo_fast, mul_mod,
    },
};

/// perform ECDSA signing on message `P` % `r` over secret key `secret_key` % `r` and nonce `k` % `r`
/// with prime subgroup generator `x, y` % `q`
pub fn ecdsa_sign<
    const NB: usize,
    P: DecomposableInto<u64> + RecomposableFrom<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    sk: &RadixCiphertext,
    k: &RadixCiphertext,
    message: P,
    generator: (P, P),
    q_modulo: P,
    r_modulo: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext) {
    // (x, y) = k * G
    println!("Calculating (x, y) = k * G");
    let (x_proj, y_proj, z_proj) = group_projective_scalar_mul::<NB, _>(
        &server_key.create_trivial_radix(generator.0, NB),
        &server_key.create_trivial_radix(generator.1, NB),
        &server_key.create_trivial_radix(1, NB),
        &k,
        q_modulo,
        server_key,
    );
    let (x, y) =
        group_projective_into_affine::<NB, _>(&x_proj, &y_proj, &z_proj, q_modulo, server_key);
    read_client_key(|client_key| {
        println!("x = {}", format(client_key.decrypt_radix::<P>(&x)));
        println!("y = {}", format(client_key.decrypt_radix::<P>(&y)));
    });
    // r = x
    // s = k^-1 * (m + r * sk)
    let r = modulo_fast::<NB, _>(&x, r_modulo, server_key);
    let k_inv = inverse_mod::<NB, _>(&k, r_modulo, server_key);
    read_client_key(|client_key| {
        println!("k^-1 = {}", format(client_key.decrypt_radix::<P>(&k_inv)));
    });
    let mrsk = add_mod::<NB, _>(
        &server_key.create_trivial_radix(message, NB),
        &mul_mod::<NB, _>(&r, &sk, r_modulo, server_key),
        r_modulo,
        server_key,
    );
    let s = mul_mod::<NB, _>(&k_inv, &mrsk, r_modulo, server_key);
    read_client_key(|client_key| {
        println!("r = {}", format(client_key.decrypt_radix::<P>(&r)));
        println!("s = {}", format(client_key.decrypt_radix::<P>(&s)));
    });

    (r, s)
}

pub fn ecdsa_verify<
    const NB: usize,
    P: DecomposableInto<u64> + RecomposableFrom<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    public_key: (RadixCiphertext, RadixCiphertext),
    signature: (RadixCiphertext, RadixCiphertext),
    message: P,
    generator: (P, P),
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    // s^-1
    let s_inv = inverse_mod::<NB, _>(&signature.1, p, server_key);
    // u1 = m * s^-1
    let u1 = mul_mod::<NB, _>(
        &s_inv,
        &server_key.create_trivial_radix(message, NB),
        p,
        server_key,
    );
    // u2 = r * s^-1
    let u2 = mul_mod::<NB, _>(&s_inv, &signature.0, p, server_key);
    // (x, y) = u1 * G + u2 * Q
    let (x_proj_1, y_proj_1, z_proj_1) = group_projective_scalar_mul::<NB, _>(
        &server_key.create_trivial_radix(generator.0, NB),
        &server_key.create_trivial_radix(generator.0, NB),
        &server_key.create_trivial_radix(1, NB),
        &u1,
        p,
        server_key,
    );
    let (x_proj_2, y_proj_2, z_proj_2) = group_projective_scalar_mul::<NB, _>(
        &public_key.0,
        &public_key.1,
        &server_key.create_trivial_radix(1, NB),
        &u2,
        p,
        server_key,
    );
    let (x_proj, y_proj, mut z_proj) = group_projective_add_projective::<NB, _>(
        &x_proj_1, &y_proj_1, &z_proj_1, &x_proj_2, &y_proj_2, &z_proj_2, p, server_key,
    );
    let mut is_z_zero = server_key.smart_scalar_eq_parallelized(&mut z_proj, 0);
    let (mut x, _y) =
        group_projective_into_affine::<NB, _>(&x_proj, &y_proj, &z_proj, p, server_key);
    let mut is_x_eq_r = server_key.smart_eq_parallelized(&mut x, &mut signature.0.clone());

    // valid if z != 0 && x == r
    server_key.smart_bitand_parallelized(&mut is_z_zero, &mut is_x_eq_r)
}
