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
    ops::{add_mod, mul_mod, mul_mod_naive},
};

pub mod helper;
pub mod ops;

fn main() {
    let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);

    //const NUM_BLOCK: usize = 8;
    //type Integer = u16;
    //let p: Integer = 13841;
    //let msg1: Integer = 13820;

    const NUM_BLOCK: usize = 128;
    type Integer = U256;
    let p = {
        let mut p = U256::ZERO;
        p.copy_from_be_byte_slice(
            &hex::decode("fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f")
                .unwrap(),
        );
        p
    };
    let msg1 = {
        let mut msg = U256::ZERO;
        //msg.copy_from_le_byte_slice(&thread_rng().gen::<[u8; 32]>());
        msg.copy_from_le_byte_slice(&[
            240, 212, 41, 181, 222, 225, 129, 5, 46, 208, 244, 180, 180, 209, 208, 121, 113, 68,
            162, 213, 221, 208, 27, 54, 57, 126, 74, 79, 45, 245, 192, 151,
        ]);
        msg % p
    };

    let ct_1 = client_key.encrypt_radix(msg1, NUM_BLOCK);
    assert_eq!(msg1, client_key.decrypt_radix::<Integer>(&ct_1));

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

    let now = Instant::now();
    let res = mul_mod_naive::<NUM_BLOCK, _>(&ct_1, &ct_1, p, &server_key);
    let elasped = now.elapsed();
    let res = client_key.decrypt_radix::<Integer>(&res);
    print!(
        "{} * {} mod {} -> {}",
        format(msg1),
        format(msg1),
        format(p),
        format(res)
    );
    println!("in {} s\n", elasped.as_secs());
}
