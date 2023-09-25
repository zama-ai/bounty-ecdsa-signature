use fhe::{
    ecdsa::{ecdsa_sign, ecdsa_sign_native, ecdsa_verify_native},
    helper::{set_client_key, u256_from_decimal_string},
    numeral::Numeral,
    ops::{
        group_jacobian::{group_projective_into_affine_native, group_projective_scalar_mul_native},
        secp256k1::prelude::*,
    },
};
use tfhe::{
    integer::{keycache::IntegerKeyCache, U256},
    shortint::prelude::PARAM_MESSAGE_2_CARRY_2,
};

fn main() {
    const NUM_BLOCK: usize = 128;
    let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);
    set_client_key(&client_key);

    let sk = u256_from_decimal_string(
        "32670510020758816978083085130507043184471273380659243275938904335757337482424",
    );
    let nonce = u256_from_decimal_string(
        "158972629851468960855479098042189567798917817837573660423710583832714848",
    );
    let msg = u256_from_decimal_string(
        "65108744961846543415519418389643270459525907322081164366671650776835723265410",
    );

    let signature_native = ecdsa_sign_native(sk, nonce, msg, *GENERATOR, *FQ_MODULO, *FR_MODULO);
    println!(
        "Native signature r: {}, s: {}",
        signature_native.0.format(),
        signature_native.1.format()
    );

    let sk_enc = client_key.encrypt_radix(sk, NUM_BLOCK);
    let nonce_enc = client_key.encrypt_radix(nonce, NUM_BLOCK);

    let signature = ecdsa_sign::<NUM_BLOCK, _>(
        &sk_enc,
        &nonce_enc,
        msg,
        *GENERATOR,
        *FQ_MODULO,
        *FR_MODULO,
        &server_key,
    );

    let signature = (
        U256::decrypt(&signature.0, &client_key),
        U256::decrypt(&signature.1, &client_key),
    );
    println!(
        "Signature: r: {}, s: {}",
        signature.0.format(),
        signature.1.format()
    );

    let public_key_projective =
        group_projective_scalar_mul_native(GENERATOR.0, GENERATOR.1, sk, *FQ_MODULO);
    let public_key = group_projective_into_affine_native(
        public_key_projective.0,
        public_key_projective.1,
        public_key_projective.2,
        *FQ_MODULO,
    );
    let is_valid = ecdsa_verify_native(
        signature, msg, public_key, *GENERATOR, *FQ_MODULO, *FR_MODULO,
    );
    println!("Is signature valid?: {}", is_valid);
}
