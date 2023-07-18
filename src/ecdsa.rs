use std::time::Instant;

use rand::Rng;
use tfhe::integer::{
    block_decomposition::{DecomposableInto, RecomposableFrom},
    RadixCiphertext, ServerKey,
};

use crate::{
    helper::{format, read_client_key},
    numeral::Numeral,
    ops::{
        add_mod,
        group_jacobian::{
            group_projective_add_projective, group_projective_into_affine,
            group_projective_scalar_mul_constant, group_projective_scalar_mul_constant_windowed,
        },
        inverse_mod, modulo_fast, mul_mod,
    },
};

/// perform ECDSA signing on message `P` % `r` over secret key `secret_key` % `r` and nonce `k` % `r`
/// with prime subgroup generator `x, y` % `q`
pub fn ecdsa_sign<const NB: usize, P: Numeral>(
    sk: &RadixCiphertext,
    k: &RadixCiphertext,
    message: P,
    generator: (P, P),
    q_modulo: P,
    r_modulo: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext) {
    // (x, y) = k * G
    #[cfg(feature = "high_level_timing")]
    let ops_start = Instant::now();
    #[cfg(feature = "high_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);
    #[cfg(feature = "high_level_timing")]
    println!("ecdsa sign start -- ref {}", task_ref);

    println!("Calculating (x, y) = k * G");
    let (x_proj, y_proj, z_proj) = group_projective_scalar_mul_constant_windowed::<8, NB, _>(
        generator.0,
        generator.1,
        k,
        q_modulo,
        server_key,
    );
    let (x, y) =
        group_projective_into_affine::<NB, _>(&x_proj, &y_proj, &z_proj, q_modulo, server_key);
    read_client_key(|client_key| {
        println!("x = {}", P::decrypt(&x, client_key).format());
        println!("y = {}", P::decrypt(&y, client_key).format());
    });
    // r = x
    // s = k^-1 * (m + r * sk)
    let r = modulo_fast::<NB, _>(&x, r_modulo, server_key);
    let k_inv = inverse_mod::<NB, _>(k, r_modulo, server_key);
    read_client_key(|client_key| {
        println!("k^-1 = {}", P::decrypt(&k_inv, client_key).format());
    });
    let mrsk = add_mod::<NB, _>(
        &server_key.create_trivial_radix(message, NB),
        &mul_mod::<NB, _>(&r, sk, r_modulo, server_key),
        r_modulo,
        server_key,
    );
    let s = mul_mod::<NB, _>(&k_inv, &mrsk, r_modulo, server_key);
    read_client_key(|client_key| {
        println!("r = {}", P::decrypt(&r, client_key).format());
        println!("s = {}", P::decrypt(&s, client_key).format());
    });

    #[cfg(feature = "high_level_timing")]
    println!(
        "ecdsa sign end, done in {:.2}s -- ref {}",
        ops_start.elapsed().as_secs_f64(),
        task_ref
    );

    (r, s)
}

pub fn ecdsa_verify<const NB: usize, P: Numeral>(
    public_key: (P, P),
    mut signature: (RadixCiphertext, RadixCiphertext),
    message: P,
    generator: (P, P),
    q_modulo: P,
    r_modulo: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    #[cfg(feature = "high_level_timing")]
    let ops_start = Instant::now();
    #[cfg(feature = "high_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);
    #[cfg(feature = "high_level_timing")]
    println!("ecdsa verify start -- ref {}", task_ref);
    // s^-1
    let s_inv = inverse_mod::<NB, _>(&signature.1, r_modulo, server_key);
    // u1 = m * s^-1
    let u1 = mul_mod::<NB, _>(
        &s_inv,
        &server_key.create_trivial_radix(message, NB),
        r_modulo,
        server_key,
    );
    // u2 = r * s^-1
    let u2 = mul_mod::<NB, _>(&s_inv, &signature.0, r_modulo, server_key);
    // (x, y) = u1 * G + u2 * Q
    let (x_proj_1, y_proj_1, z_proj_1) = group_projective_scalar_mul_constant_windowed::<8, NB, _>(
        generator.0,
        generator.1,
        &u1,
        q_modulo,
        server_key,
    );
    let (x_proj_2, y_proj_2, z_proj_2) = group_projective_scalar_mul_constant_windowed::<8, NB, _>(
        public_key.0,
        public_key.1,
        &u2,
        q_modulo,
        server_key,
    );
    let (x_proj, y_proj, mut z_proj) = group_projective_add_projective::<NB, _>(
        &x_proj_1, &y_proj_1, &z_proj_1, &x_proj_2, &y_proj_2, &z_proj_2, q_modulo, server_key,
    );
    let mut is_z_zero = server_key.smart_scalar_eq_parallelized(&mut z_proj, 0);
    let (x, _y) =
        group_projective_into_affine::<NB, _>(&x_proj, &y_proj, &z_proj, q_modulo, server_key);
    let mut x_mod_scalar = modulo_fast::<NB, _>(&x, r_modulo, server_key);
    let mut is_x_eq_r = server_key.smart_eq_parallelized(&mut x_mod_scalar, &mut signature.0);

    #[cfg(feature = "high_level_timing")]
    println!(
        "ecdsa verify done in {:.2}s -- ref {}",
        ops_start.elapsed().as_secs_f64(),
        task_ref
    );
    // valid if z != 0 && x == r
    server_key.smart_bitand_parallelized(&mut is_z_zero, &mut is_x_eq_r)
}
