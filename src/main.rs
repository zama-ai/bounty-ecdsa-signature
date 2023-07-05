#![allow(unused_imports)]
#![feature(iter_array_chunks)]

use std::time::Instant;

use num_bigint::BigInt;
use rand::{thread_rng, Rng};
use tfhe::{
    integer::{keycache::IntegerKeyCache, U256},
    shortint::parameters::PARAM_MESSAGE_2_CARRY_2,
};

use crate::{
    ecdsa::ecdsa_sign,
    helper::{format, u256_from_decimal_string},
    ops::{
        add_mod, double_mod,
        group::{
            group_projective_add_projective, group_projective_double, group_projective_into_affine,
            group_projective_scalar_mul,
        },
        inverse_mod, mul_mod, mul_mod_bitwise, pow_mod,
    },
};

pub mod ecdsa;
pub mod helper;
pub mod ops;

fn main() {
    let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);

    //const NUM_BLOCK: usize = 128;
    //type Integer = U256;
    //let q: Integer = u256_from_decimal_string(
    //"115792089237316195423570985008687907853269984665640564039457584007908834671663",
    //);
    //let r: Integer = u256_from_decimal_string(
    //"115792089237316195423570985008687907852837564279074904382605163141518161494337",
    //);
    //let g_x: Integer = u256_from_decimal_string(
    //"55066263022277343669578718895168534326250603453777594175500187360389116729240",
    //);
    //let g_y: Integer = u256_from_decimal_string(
    //"32670510020758816978083085130507043184471273380659243275938904335757337482424",
    //);
    //let secret_key: Integer = u256_from_decimal_string(
    //"50374736940266874356472946887032805385643216083156316185256796863771578718026",
    //);
    //let nonce: Integer = u256_from_decimal_string(
    //"31534706826920421062170617943943122326357919410239174025491525139969337821575",
    //);
    //let message: Integer = u256_from_decimal_string(
    //"85842141585962005810830452070980593281629797958423199730619494346446023220882",
    //);

    //let generator = (g_x, g_y);
    //let secret_key_enc = client_key.encrypt_radix(secret_key, NUM_BLOCK);
    //let nonce_enc = client_key.encrypt_radix(nonce, NUM_BLOCK);
    //let now = Instant::now();
    //println!("Start ecdsa signing");
    //let (r, s) = ecdsa_sign::<NUM_BLOCK, _>(
    //&secret_key_enc,
    //&nonce_enc,
    //message,
    //generator,
    //q,
    //r,
    //&server_key,
    //&client_key,
    //);
    //let elasped = now.elapsed();
    //println!(
    //"Signed message {} with {} and got signature r = {}, s = {}",
    //format(message),
    //format(secret_key),
    //format(client_key.decrypt_radix::<Integer>(&r)),
    //format(client_key.decrypt_radix::<Integer>(&s))
    //);
    //println!("Signed in {:.2} s\n", elasped.as_secs_f32());

    const NUM_BLOCK: usize = 4;
    type Integer = u8;
    let p: Integer = 251;
    let x1: Integer = 8;
    let y1: Integer = 45;
    let x2: Integer = 26;
    let y2: Integer = 55;

    //const NUM_BLOCK: usize = 8;
    //type Integer = u16;
    //let p: Integer = 65521;
    //let x: Integer = 13897;
    //let y: Integer = 62399;

    //const NUM_BLOCK: usize = 128;
    //type Integer = U256;
    //let p = {
    //let mut p = U256::ZERO;
    //p.copy_from_be_byte_slice(
    //&hex::decode("fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f")
    //.unwrap(),
    //);
    //p
    //};
    //let msg1 = {
    //let mut msg = U256::ZERO;
    ////msg.copy_from_le_byte_slice(&thread_rng().gen::<[u8; 32]>());
    //msg.copy_from_le_byte_slice(&[
    //240, 212, 41, 181, 222, 225, 129, 5, 46, 208, 244, 180, 180, 209, 208, 121, 113, 68,
    //162, 213, 221, 208, 27, 54, 57, 126, 74, 79, 45, 245, 192, 151,
    //]);
    //msg % p
    //};

    //let add_mod_int = |x: Integer, y: Integer| -> Integer {
    //((BigInt::from(x) + BigInt::from(y)) % &BigInt::from(p))
    //.try_into()
    //.unwrap()
    //};
    //let mul_mod_int = |x: Integer, y: Integer| -> Integer {
    //((BigInt::from(x) * BigInt::from(y)) % &BigInt::from(p))
    //.try_into()
    //.unwrap()
    //};

    let ct_x1 = client_key.encrypt_radix(x1, NUM_BLOCK);
    let ct_y1 = client_key.encrypt_radix(y1, NUM_BLOCK);
    let ct_x2 = client_key.encrypt_radix(x2, NUM_BLOCK);
    let ct_y2 = client_key.encrypt_radix(y2, NUM_BLOCK);
    assert_eq!(x1, client_key.decrypt_radix::<Integer>(&ct_x1));
    assert_eq!(y1, client_key.decrypt_radix::<Integer>(&ct_y1));
    assert_eq!(x2, client_key.decrypt_radix::<Integer>(&ct_x2));
    assert_eq!(y2, client_key.decrypt_radix::<Integer>(&ct_y2));

    //let now = Instant::now();
    //let res = add_mod::<NUM_BLOCK, _>(&ct_x1, &ct_y1, p, &server_key);
    //let elasped = now.elapsed();
    //let res = client_key.decrypt_radix::<Integer>(&res);
    //println!(
    //"{} + {} mod {} -> {}",
    //format(x1),
    //format(y1),
    //format(p),
    //format(res)
    //);
    //println!("should be {}", add_mod_int(x1, y1));
    //println!("add mod in {:.2} s\n", elasped.as_secs_f32());

    //let now = Instant::now();
    //let res = double_mod::<NUM_BLOCK, _>(&ct_x1, p, &server_key);
    //let elasped = now.elapsed();
    //let res = client_key.decrypt_radix::<Integer>(&res);
    //println!("{} * 2 mod {} -> {}", format(x1), format(p), format(res));
    //println!("should be {}", add_mod_int(x1, x1));
    //println!("add mod in {:.2} s\n", elasped.as_secs_f32());

    //let now = Instant::now();
    //let res = mul_mod::<NUM_BLOCK, _>(&ct_x1, &ct_y1, p, &server_key);
    //let elasped = now.elapsed();
    //let res = client_key.decrypt_radix::<Integer>(&res);
    //println!(
    //"{} * {} mod {} -> {}",
    //format(x1),
    //format(y1),
    //format(p),
    //format(res)
    //);
    //println!("should be {}", mul_mod_int(x1, y1));
    //println!("mul mod in {:.2} s\n", elasped.as_secs_f32());

    //let now = Instant::now();
    //let res = pow_mod::<NUM_BLOCK, _>(&ct_x1, &ct_y1, p, &server_key);
    //let elasped = now.elapsed();
    //println!(
    //"{}^{} % {} -> {}",
    //format(x1),
    //format(y1),
    //format(p),
    //format(client_key.decrypt_radix::<Integer>(&res)),
    //);
    //println!("pow mod in {:.2} s\n", elasped.as_secs_f32());

    //let now = Instant::now();
    //let res = inverse_mod::<NUM_BLOCK, _>(&ct_x1, p, &server_key);
    //let res_decoded = client_key.decrypt_radix::<Integer>(&res);
    //println!(
    //"{}^-1 % {} -> {}",
    //format(x1),
    //format(p),
    //format(res_decoded)
    //);
    //let elasped = now.elapsed();
    //println!("inverse mod in {:.2} s\n", elasped.as_secs_f32());

    //let now = Instant::now();
    //let (x_new, y_new, z_new) = group_projective_double::<NUM_BLOCK, _>(
    //&ct_x1,
    //&ct_y1,
    //&server_key.create_trivial_radix(1, NUM_BLOCK),
    //p,
    //&server_key,
    //);
    //let elasped = now.elapsed();
    //let x_dec = client_key.decrypt_radix::<Integer>(&x_new);
    //let y_dec = client_key.decrypt_radix::<Integer>(&y_new);
    //let z_dec = client_key.decrypt_radix::<Integer>(&z_new);
    //println!(
    //"{},{},{} * 2 -> {},{},{}",
    //format(x1),
    //format(y1),
    //format(1),
    //format(x_dec),
    //format(y_dec),
    //format(z_dec)
    //);
    //println!("group double in {:.2} s", elasped.as_secs_f32());

    let now = Instant::now();
    let (x_new, y_new, z_new) = group_projective_add_projective::<NUM_BLOCK, _>(
        &ct_x1,
        &ct_y1,
        &client_key.encrypt_radix(1, NUM_BLOCK),
        &ct_x2,
        &ct_y2,
        &client_key.encrypt_radix(1, NUM_BLOCK),
        p,
        &server_key,
    );
    let elasped = now.elapsed();
    let x_dec = client_key.decrypt_radix::<Integer>(&x_new);
    let y_dec = client_key.decrypt_radix::<Integer>(&y_new);
    let z_dec = client_key.decrypt_radix::<Integer>(&z_new);
    println!(
        "{},{},{} + {},{},{} -> {},{},{}",
        format(x1),
        format(y1),
        format(1),
        format(x2),
        format(y2),
        format(1),
        format(x_dec),
        format(y_dec),
        format(z_dec)
    );
    println!("group add in {} s", elasped.as_secs_f32());

    //let now = Instant::now();
    //let (x_new, y_new, z_new) = group_projective_scalar_mul::<NUM_BLOCK, _>(
    //&ct_x1,
    //&ct_y1,
    //&client_key.encrypt_radix(1, NUM_BLOCK),
    //&ct_x2,
    //p,
    //&server_key,
    //&client_key,
    //);
    //let x_dec = client_key.decrypt_radix::<Integer>(&x_new);
    //let y_dec = client_key.decrypt_radix::<Integer>(&y_new);
    //let z_dec = client_key.decrypt_radix::<Integer>(&z_new);
    //let elasped = now.elapsed();
    //println!(
    //"{},{},{} * {} -> {},{},{}",
    //format(x1),
    //format(y1),
    //format(1),
    //format(x2),
    //format(x_dec),
    //format(y_dec),
    //format(z_dec)
    //);
    //println!("group scalar mul in {:.2} s", elasped.as_secs_f32());

    let now = Instant::now();
    let (x_aff, y_aff) =
        group_projective_into_affine::<NUM_BLOCK, _>(&x_new, &y_new, &z_new, p, &server_key);
    let elasped = now.elapsed();
    println!(
        "{},{},{} -> {},{}",
        format(x_dec),
        format(y_dec),
        format(z_dec),
        format(client_key.decrypt_radix::<Integer>(&x_aff)),
        format(client_key.decrypt_radix::<Integer>(&y_aff))
    );
    println!(
        "group projective into affine in {:.2}s",
        elasped.as_secs_f32()
    );
}
