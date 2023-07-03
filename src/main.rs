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
    ops::{add_mod, group_projective_double, mul_mod, pow_mod},
};

pub mod helper;
pub mod ops;

fn main() {
    let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);

    const NUM_BLOCK: usize = 4;
    type Integer = u8;
    let p: Integer = 251;
    let x: Integer = 18;
    let y: Integer = 130;

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
    //let res = add_mod::<NUM_BLOCK, _>(&ct_1, &ct_1, p, &server_key);
    //let elasped = now.elapsed();
    //let res = client_key.decrypt_radix::<Integer>(&res);
    //println!(
    //"{} + {} mod {} -> {}",
    //format(msg1),
    //format(msg1),
    //format(p),
    //format(res)
    //);
    //println!("in {} s\n", elasped.as_secs());

    //let now = Instant::now();
    //let res = mul_mod::<NUM_BLOCK, _>(&ct_1, &ct_1, p, &server_key);
    //let elasped = now.elapsed();
    //let res = client_key.decrypt_radix::<Integer>(&res);
    //print!(
    //"{} * {} mod {} -> {}",
    //format(msg1),
    //format(msg1),
    //format(p),
    //format(res)
    //);
    //println!("in {} s\n", elasped.as_secs());

    //let now = Instant::now();
    //let (x_new, y_new, z_new) =
    //group_projective_double::<NUM_BLOCK, _>(&ct_x, &ct_y, &ct_z, p, &server_key);
    //let elasped = now.elapsed();
    //print!(
    //"{},{},{} * 2 -> {},{},{}",
    //format(x),
    //format(y),
    //format(z),
    //format(client_key.decrypt_radix::<Integer>(&x_new)),
    //format(client_key.decrypt_radix::<Integer>(&y_new)),
    //format(client_key.decrypt_radix::<Integer>(&z_new))
    //);
    //println!("in {} s\n", elasped.as_secs());

    let now = Instant::now();
    let res = pow_mod::<NUM_BLOCK, _>(&ct_x, &ct_y, p, &server_key);
    let elasped = now.elapsed();
    print!(
        "{}^{} % {} -> {}",
        format(x),
        format(y),
        format(p),
        format(client_key.decrypt_radix::<Integer>(&res)),
    );
    println!("in {} s\n", elasped.as_secs());
}
