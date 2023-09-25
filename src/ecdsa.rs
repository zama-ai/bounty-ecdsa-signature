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
            group_projective_add_projective, group_projective_add_projective_native,
            group_projective_into_affine, group_projective_into_affine_inv,
            group_projective_into_affine_native, group_projective_scalar_mul_constant,
            group_projective_scalar_mul_constant_windowed, group_projective_scalar_mul_native,
        },
        inverse_mod, inverse_mods,
        mersenne::mod_mersenne,
        modulo_fast, mul_mod,
        native::{add_mod_native, inverse_mod_native, modulo_native, mul_mod_native},
    },
    WINDOW,
};

/// perform homomorphic ECDSA signing on message `P` % `r` over secret key `secret_key` % `r` and nonce `k` % `r`
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
    println!("ECDSA sign start");
    println!("Calculating (x, y) = k * G");
    let ops_start = Instant::now();
    let (x_proj, y_proj, z_proj) = group_projective_scalar_mul_constant_windowed::<WINDOW, NB, _>(
        generator.0,
        generator.1,
        k,
        q_modulo,
        server_key,
    );
    let (z_inv, k_inv) = rayon::join(
        || inverse_mod::<NB, _>(&z_proj, q_modulo, server_key),
        || inverse_mod::<NB, _>(k, r_modulo, server_key),
    );
    let (x, y) =
        group_projective_into_affine_inv::<NB, _>(&x_proj, &y_proj, &z_inv, q_modulo, server_key);
    read_client_key(|client_key| {
        println!("x = {}", P::decrypt(&x, client_key).format());
        println!("y = {}", P::decrypt(&y, client_key).format());
    });
    // r = x
    // s = k^-1 * (m + r * sk)
    let r = if q_modulo > r_modulo && q_modulo <= P::TWO * r_modulo {
        modulo_fast::<NB, _>(&x, r_modulo, server_key)
    } else {
        mod_mersenne::<NB, _>(&x, r_modulo, server_key)
    };
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

    println!(
        "ECDSA sign end, done in {:.2}s",
        ops_start.elapsed().as_secs_f64(),
    );

    (r, s)
}

/// verify ECDSA signature 
pub fn ecdsa_sign_native<P: Numeral>(
    sk: P,
    k: P,
    message: P,
    generator: (P, P),
    q_modulo: P,
    r_modulo: P,
) -> (P, P) {
    // (x, y) = k * G
    let (x, y, z) = group_projective_scalar_mul_native(generator.0, generator.1, k, q_modulo);
    let (x, _y) = group_projective_into_affine_native(x, y, z, q_modulo);
    // r = x
    // s = k^-1 * (m + r * sk)
    let r = modulo_native(x, r_modulo);
    let k_inv = inverse_mod_native(k, r_modulo);
    let mrsk = add_mod_native(message, mul_mod_native(r, sk, r_modulo), r_modulo);
    let s = mul_mod_native(k_inv, mrsk, r_modulo);

    (r, s)
}

pub fn ecdsa_verify_native<P: Numeral>(
    signature: (P, P),
    message: P,
    public_key: (P, P),
    generator: (P, P),
    q_modulo: P,
    r_modulo: P,
) -> bool {
    if signature.0 >= r_modulo || signature.1 >= r_modulo {
        return false;
    }

    let s_inv = inverse_mod_native(signature.1, r_modulo);
    let u1 = mul_mod_native(message, s_inv, r_modulo);
    let u2 = mul_mod_native(signature.0, s_inv, r_modulo);
    let (x1, y1, z1) = group_projective_scalar_mul_native(generator.0, generator.1, u1, q_modulo);
    let (x2, y2, z2) = group_projective_scalar_mul_native(public_key.0, public_key.1, u2, q_modulo);
    let (x, y, z) = group_projective_add_projective_native(x1, y1, z1, x2, y2, z2, q_modulo);
    let (x, _y) = group_projective_into_affine_native(x, y, z, q_modulo);

    return signature.0 == modulo_native(x, r_modulo);
}

#[cfg(test)]
mod tests {
    use tfhe::{integer::keycache::IntegerKeyCache, shortint::prelude::PARAM_MESSAGE_2_CARRY_2};

    use crate::{
        ecdsa::ecdsa_sign,
        helper::set_client_key,
        numeral::Numeral,
        ops::group_jacobian::{
            group_projective_double_native, group_projective_into_affine_native,
            group_projective_scalar_mul_native,
        },
    };

    use super::{ecdsa_sign_native, ecdsa_verify_native};

    #[test]
    fn correct_ecdsa_sign_verify_native() {
        let q_modulo: u8 = 211;
        let gx: u8 = 4;
        let gy: u8 = 156;
        let r_modulo: u8 = 199;

        let sk = 111;
        let k = 71;
        let message = 89;
        let pk_projective = group_projective_scalar_mul_native(gx, gy, sk, q_modulo);
        let pk = group_projective_into_affine_native(
            pk_projective.0,
            pk_projective.1,
            pk_projective.2,
            q_modulo,
        );
        let (r, s) = ecdsa_sign_native(sk, k, message, (gx, gy), q_modulo, r_modulo);
        let is_valid = ecdsa_verify_native((r, s), message, pk, (gx, gy), q_modulo, r_modulo);
        assert!(is_valid, "ECDSA signature is invalid");
    }

    #[test]
    fn correct_ecdsa_sign_verify() {
        let q_modulo: u8 = 211;
        let gx: u8 = 4;
        let gy: u8 = 156;
        let r_modulo: u8 = 199;

        let sk = 111;
        let k = 71;
        let message = 89;
        let (rx, ry) = ecdsa_sign_native(sk, k, message, (gx, gy), q_modulo, r_modulo);

        let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);
        const NUM_BLOCK: usize = 4;

        let enc_sk = client_key.encrypt_radix(sk, NUM_BLOCK);
        let enc_k = client_key.encrypt_radix(k, NUM_BLOCK);

        let (enc_rx, enc_ry) = ecdsa_sign::<NUM_BLOCK, _>(
            &enc_sk,
            &enc_k,
            message,
            (gx, gy),
            q_modulo,
            r_modulo,
            &server_key,
        );

        assert_eq!(rx, u8::decrypt(&enc_rx, &client_key));
        assert_eq!(ry, u8::decrypt(&enc_ry, &client_key));

        let pk_projective = group_projective_scalar_mul_native(gx, gy, sk, q_modulo);
        let pk = group_projective_into_affine_native(
            pk_projective.0,
            pk_projective.1,
            pk_projective.2,
            q_modulo,
        );
        let is_valid = ecdsa_verify_native((rx, ry), message, pk, (gx, gy), q_modulo, r_modulo);
        assert!(is_valid, "ECDSA signature is invalid");
    }
}
