#![feature(iter_array_chunks)]
#![feature(result_option_inspect)]
#![allow(unused_imports)]

use std::{cell::Cell, sync::RwLock, time::Instant};

use lazy_static::lazy_static;
use num_bigint::BigInt;
use rand::{thread_rng, Rng};
use tfhe::{
    integer::{keycache::IntegerKeyCache, ClientKey, U256},
    shortint::parameters::PARAM_MESSAGE_2_CARRY_2,
};

#[cfg(not(feature = "jacobian"))]
use crate::ops::group_homogenous::{
    group_projective_add_projective, group_projective_double, group_projective_into_affine,
    group_projective_scalar_mul,
};
#[cfg(feature = "jacobian")]
use crate::ops::group_jacobian::{
    group_projective_add_projective, group_projective_double, group_projective_into_affine,
    group_projective_scalar_mul,
};
use crate::{
    ecdsa::ecdsa_sign,
    helper::{format, from_bigint, to_bigint, u256_from_decimal_string},
    ops::{add_mod, double_mod, inverse_mod, mul_mod, mul_mod_bitwise, pow_mod},
};

pub mod ecdsa;
pub mod field;
pub mod helper;
pub mod ops;

lazy_static! {
    pub static ref CLIENT_KEY: RwLock<Option<ClientKey>> = RwLock::new(None);
}

fn main() {
    let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);
    *CLIENT_KEY.write().unwrap() = Some(client_key.clone());

    #[cfg(feature = "go_big")]
    const NUM_BLOCK: usize = 128;
    #[cfg(feature = "go_big")]
    type Integer = U256;
    #[cfg(feature = "go_big")]
    let (p, x1, y1, x2, y2) = {
        let p: Integer = u256_from_decimal_string(
            "115792089237316195423570985008687907853269984665640564039457584007908834671663",
        );
        let x1: Integer = u256_from_decimal_string(
            "55066263022277343669578718895168534326250603453777594175500187360389116729240",
        );
        let y1: Integer = u256_from_decimal_string(
            "32670510020758816978083085130507043184471273380659243275938904335757337482424",
        );
        let x2: Integer = u256_from_decimal_string(
            "55066263022277343669578718895168534326250603453777594175500187360389116729240",
        );
        let y2: Integer = u256_from_decimal_string(
            "32670510020758816978083085130507043184471273380659243275938904335757337482424",
        );
        (p, x1, y1, x2, y2)
    };

    #[cfg(not(feature = "go_big"))]
    const NUM_BLOCK: usize = 8;
    #[cfg(not(feature = "go_big"))]
    type Integer = u16;

    #[cfg(not(feature = "go_big"))]
    let (p, x1, y1, x2, y2) = {
        let p: Integer = 65535;
        let x1: Integer = 50725;
        let y1: Integer = 64006;
        let x2: Integer = 34884;
        let y2: Integer = 48022;
        (p, x1, y1, x2, y2)
    };

    let add_mod_int = |x: Integer, y: Integer| -> Integer {
        let x = to_bigint(x);
        let y = to_bigint(y);
        let p = to_bigint(p);
        from_bigint(&((&x + &y) % &p))
    };
    let mul_mod_int = |x: Integer, y: Integer| -> Integer {
        let x = to_bigint(x);
        let y = to_bigint(y);
        let p = to_bigint(p);
        from_bigint(&((&x * &y) % &p))
    };

    let ct_x1 = client_key.encrypt_radix(x1, NUM_BLOCK);
    let ct_y1 = client_key.encrypt_radix(y1, NUM_BLOCK);
    let ct_x2 = client_key.encrypt_radix(x2, NUM_BLOCK);
    let ct_y2 = client_key.encrypt_radix(y2, NUM_BLOCK);
    assert_eq!(x1, client_key.decrypt_radix::<Integer>(&ct_x1));
    assert_eq!(y1, client_key.decrypt_radix::<Integer>(&ct_y1));
    assert_eq!(x2, client_key.decrypt_radix::<Integer>(&ct_x2));
    assert_eq!(y2, client_key.decrypt_radix::<Integer>(&ct_y2));
    println!("Finished asserting ciphertext");

    let now = Instant::now();
    let res = add_mod::<NUM_BLOCK, _>(&ct_x1, &ct_y1, p, &server_key);
    let elasped = now.elapsed();
    let res = client_key.decrypt_radix::<Integer>(&res);
    println!(
        "{} + {} mod {} -> {}",
        format(x1),
        format(y1),
        format(p),
        format(res)
    );
    println!("should be {}", format(add_mod_int(x1, y1)));
    println!("add mod in {:.2} s\n", elasped.as_secs_f32());

    let now = Instant::now();
    let res = double_mod::<NUM_BLOCK, _>(&ct_x1, p, &server_key);
    let elasped = now.elapsed();
    let res = client_key.decrypt_radix::<Integer>(&res);
    println!("{} * 2 mod {} -> {}", format(x1), format(p), format(res));
    println!("should be {}", format(add_mod_int(x1, x1)));
    println!("double mod in {:.2} s\n", elasped.as_secs_f32());

    let now = Instant::now();
    let res = mul_mod::<NUM_BLOCK, _>(&ct_x1, &ct_y1, p, &server_key);
    let elasped = now.elapsed();
    let res = client_key.decrypt_radix::<Integer>(&res);
    println!(
        "{} * {} mod {} -> {}",
        format(x1),
        format(y1),
        format(p),
        format(res)
    );
    println!("should be {}", format(mul_mod_int(x1, y1)));
    println!("mul mod in {:.2} s\n", elasped.as_secs_f32());

    let now = Instant::now();
    let res = inverse_mod::<NUM_BLOCK, _>(&ct_x1, p, &server_key);
    let res_decoded = client_key.decrypt_radix::<Integer>(&res);
    println!(
        "{}^-1 % {} -> {}",
        format(x1),
        format(p),
        format(res_decoded)
    );
    let elasped = now.elapsed();
    println!("inverse mod in {:.2} s\n", elasped.as_secs_f32());

    // let now = Instant::now();
    // let res = pow_mod::<NUM_BLOCK, _>(&ct_x1, &ct_y1, p, &server_key);
    // let elasped = now.elapsed();
    // println!(
    //     "{}^{} % {} -> {}",
    //     format(x1),
    //     format(y1),
    //     format(p),
    //     format(client_key.decrypt_radix::<Integer>(&res)),
    // );
    // println!("pow mod in {:.2} s\n", elasped.as_secs_f32());

    let now = Instant::now();
    let (x_new, y_new, z_new) = group_projective_double::<NUM_BLOCK, _>(
        &ct_x1,
        &ct_y1,
        &server_key.create_trivial_radix(1, NUM_BLOCK),
        p,
        &server_key,
    );
    let elasped = now.elapsed();
    let x_dec = client_key.decrypt_radix::<Integer>(&x_new);
    let y_dec = client_key.decrypt_radix::<Integer>(&y_new);
    let z_dec = client_key.decrypt_radix::<Integer>(&z_new);
    println!(
        "{},{},{} * 2 -> {},{},{}",
        format(x1),
        format(y1),
        format(1),
        format(x_dec),
        format(y_dec),
        format(z_dec)
    );
    println!("group double in {:.2} s", elasped.as_secs_f32());

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

    let now = Instant::now();
    let (x_new, y_new, z_new) = group_projective_scalar_mul::<NUM_BLOCK, _>(
        &ct_x1,
        &ct_y1,
        &client_key.encrypt_radix(1, NUM_BLOCK),
        &ct_x2,
        p,
        &server_key,
    );
    let x_dec = client_key.decrypt_radix::<Integer>(&x_new);
    let y_dec = client_key.decrypt_radix::<Integer>(&y_new);
    let z_dec = client_key.decrypt_radix::<Integer>(&z_new);
    let elasped = now.elapsed();
    println!(
        "{},{},{} * {} -> {},{},{}",
        format(x1),
        format(y1),
        format(1),
        format(x2),
        format(x_dec),
        format(y_dec),
        format(z_dec)
    );
    println!("group scalar mul in {:.2} s", elasped.as_secs_f32());

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
