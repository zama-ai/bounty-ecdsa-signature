#![allow(unused_imports)]
#![feature(iter_array_chunks)]

use std::time::Instant;

use rand::{thread_rng, Rng};
use tfhe::{
    integer::{keycache::IntegerKeyCache, U256},
    shortint::parameters::PARAM_MESSAGE_2_CARRY_2,
};

use crate::{
    helper::format,
    ops::{
        add_mod,
        group::{
            group_projective_add_projective, group_projective_double, group_projective_into_affine,
            group_projective_scalar_mul,
        },
        inverse_mod, mul_mod, pow_mod,
    },
};

pub mod ecdsa;
pub mod helper;
pub mod ops;

fn main() {
    let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);

    const NUM_BLOCK: usize = 4;
    type Integer = u8;
    let p: Integer = 251;
    let x: Integer = 8;
    let y: Integer = 45;

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

    let ct_x = client_key.encrypt_radix(x, NUM_BLOCK);
    let ct_y = client_key.encrypt_radix(y, NUM_BLOCK);
    //let ct_z = client_key.encrypt_radix(z, NUM_BLOCK);
    assert_eq!(x, client_key.decrypt_radix::<Integer>(&ct_x));
    assert_eq!(y, client_key.decrypt_radix::<Integer>(&ct_y));
    //assert_eq!(z, client_key.decrypt_radix::<Integer>(&ct_z));

    //let now = Instant::now();
    //let res = add_mod::<NUM_BLOCK, _>(&ct_x, &ct_y, p, &server_key);
    //let elasped = now.elapsed();
    //let res = client_key.decrypt_radix::<Integer>(&res);
    //println!(
    //"{} + {} mod {} -> {}",
    //format(x),
    //format(y),
    //format(p),
    //format(res)
    //);
    //println!("in {:.2} s\n", elasped.as_secs_f32());

    //let now = Instant::now();
    //let res = mul_mod::<NUM_BLOCK, _>(&ct_x, &ct_y, p, &server_key);
    //let elasped = now.elapsed();
    //let res = client_key.decrypt_radix::<Integer>(&res);
    //println!(
    //"{} * {} mod {} -> {}",
    //format(x),
    //format(y),
    //format(p),
    //format(res)
    //);
    //println!("in {:.2} s\n", elasped.as_secs_f32());

    //let now = Instant::now();
    //let res = pow_mod::<NUM_BLOCK, _>(&ct_x, &ct_y, p, &server_key);
    //let elasped = now.elapsed();
    //println!(
    //"{}^{} % {} -> {}",
    //format(x),
    //format(y),
    //format(p),
    //format(client_key.decrypt_radix::<Integer>(&res)),
    //);
    //println!("in {:.2} s\n", elasped.as_secs_f32());

    //let now = Instant::now();
    //let res = inverse_mod::<NUM_BLOCK, _>(&ct_x, p, &server_key);
    //let res_decoded = client_key.decrypt_radix::<Integer>(&res);
    //println!(
    //"{}^-1 % {} -> {}",
    //format(x),
    //format(p),
    //format(res_decoded)
    //);
    //let elasped = now.elapsed();
    //println!("in {:.2} s\n", elasped.as_secs_f32());

    //let now = Instant::now();
    //let (x_new, y_new, z_new) = group_projective_double::<NUM_BLOCK, _>(
    //&ct_x,
    //&ct_y,
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
    //format(x),
    //format(y),
    //format(1),
    //format(x_dec),
    //format(y_dec),
    //format(z_dec)
    //);
    //println!("group double in {} s", elasped.as_secs());

    //let now = Instant::now();
    //let (x_new, y_new, z_new) = group_projective_add_projective::<NUM_BLOCK, _>(
    //&client_key.encrypt_radix(8, NUM_BLOCK),
    //&client_key.encrypt_radix(45, NUM_BLOCK),
    //&server_key.create_trivial_radix(1, NUM_BLOCK),
    //&client_key.encrypt_radix(26, NUM_BLOCK),
    //&client_key.encrypt_radix(55, NUM_BLOCK),
    //&server_key.create_trivial_radix(1, NUM_BLOCK),
    ////&server_key.create_trivial_radix(0, NUM_BLOCK),
    ////&server_key.create_trivial_radix(0, NUM_BLOCK),
    ////&server_key.create_trivial_radix(0, NUM_BLOCK),
    //p,
    //&server_key,
    //);
    //let elasped = now.elapsed();
    //let x_dec = client_key.decrypt_radix::<Integer>(&x_new);
    //let y_dec = client_key.decrypt_radix::<Integer>(&y_new);
    //let z_dec = client_key.decrypt_radix::<Integer>(&z_new);
    //println!(
    //"{},{},{} + {},{},{} -> {},{},{}",
    //format(8),
    //format(45),
    //format(1),
    //format(26),
    //format(55),
    //format(1),
    ////format(0),
    ////format(0),
    ////format(0),
    //format(x_dec),
    //format(y_dec),
    //format(z_dec)
    //);
    //println!("group double in {} s", elasped.as_secs());

    let now = Instant::now();
    let mul: Integer = 234;
    let (x_new, y_new, z_new) = group_projective_scalar_mul::<NUM_BLOCK, _>(
        &ct_x,
        &ct_y,
        &server_key.create_trivial_radix(1, NUM_BLOCK),
        &client_key.encrypt_radix(mul, NUM_BLOCK),
        p,
        &server_key,
        &client_key,
    );
    let x_dec = client_key.decrypt_radix::<Integer>(&x_new);
    let y_dec = client_key.decrypt_radix::<Integer>(&y_new);
    let z_dec = client_key.decrypt_radix::<Integer>(&z_new);
    let elasped = now.elapsed();
    println!(
        "{},{},{} * {} -> {},{},{}",
        format(x),
        format(y),
        format(1),
        format(mul),
        format(x_dec),
        format(y_dec),
        format(z_dec)
    );
    println!("group scalar mul in {:.2} s", elasped.as_secs_f64());

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
    println!("group projective into affine in {}s", elasped.as_secs());
}
