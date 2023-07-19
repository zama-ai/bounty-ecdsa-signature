#![allow(unused_imports)]

use fhe::{
    helper::{format, set_client_key, u256_from_decimal_string},
    numeral::Numeral,
    ops::{
        add_mod, double_mod,
        group_jacobian::{
            group_projective_into_affine, group_projective_scalar_mul_constant_windowed,
        },
        inverse_mod, inverse_mods, mul_mod,
        native::{add_mod_native, inverse_mod_native, mul_mod_native, sub_mod_native},
        sub_mod,
    },
    stats::ProtocolStats,
};
use std::time::Instant;
use tfhe::{
    integer::{keycache::IntegerKeyCache, U256},
    shortint::parameters::PARAM_MESSAGE_2_CARRY_2,
};

fn main() {
    let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);
    set_client_key(&client_key);

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

    //#[cfg(not(feature = "go_big"))]
    //const NUM_BLOCK: usize = 8;
    //#[cfg(not(feature = "go_big"))]
    //type Integer = u16;
    //#[cfg(not(feature = "go_big"))]
    //let (p, x1, y1, x2, y2) = {
    //let p: Integer = 65521;
    //let x1: Integer = 50725;
    //let y1: Integer = 64006;
    //let x2: Integer = 34884;
    //let y2: Integer = 48022;
    //(p, x1, y1, x2, y2)
    //};

    #[cfg(not(feature = "go_big"))]
    const NUM_BLOCK: usize = 4;
    #[cfg(not(feature = "go_big"))]
    type Integer = u8;
    #[cfg(not(feature = "go_big"))]
    let (p, x1, y1, x2, y2) = {
        let p: Integer = 251;
        let x1: Integer = 8;
        let y1: Integer = 45;
        let x2: Integer = 26;
        let y2: Integer = 55;
        (p, x1, y1, x2, y2)
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

    //let now = Instant::now();
    //let res = double_mod::<NUM_BLOCK, _>(&ct_x1, p, &server_key);
    //let elasped = now.elapsed();
    //let res = client_key.decrypt_radix::<Integer>(&res);
    //println!("{} * 2 mod {} -> {}", format(x1), format(p), format(res));
    //println!("should be {}", format(add_mod_native(x1, x1, p)));
    //println!("double mod in {:.2} s\n", elasped.as_secs_f32());

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
    //println!("should be {}", format(add_mod_native(x1, y1, p)));
    //println!("add mod in {:.2} s\n", elasped.as_secs_f32());

    //let now = Instant::now();
    //let res = sub_mod::<NUM_BLOCK, _>(&ct_x1, &ct_y1, p, &server_key);
    //let elasped = now.elapsed();
    //let res = client_key.decrypt_radix::<Integer>(&res);
    //println!(
    //"{} - {} mod {} -> {}",
    //format(x1),
    //format(y1),
    //format(p),
    //format(res)
    //);
    //println!("should be {}", format(sub_mod_native(x1, y1, p)));
    //println!("sub mod in {:.2} s\n", elasped.as_secs_f32());

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
    //println!("should be {}", format(mul_mod_native(x1, y1, p)));
    //println!("mul mod in {:.2} s\n", elasped.as_secs_f32());

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
    //println!("should be {}", format(inverse_mod_native(x1, p)));
    //println!("inverse mod in {:.2} s\n", elasped.as_secs_f32());

    let now = Instant::now();
    let res = inverse_mods::<NUM_BLOCK, _>(&[ct_x1.clone(), ct_y1.clone()], p, &server_key);
    let res1_decoded = Integer::decrypt(&res[0], &client_key);
    let res2_decoded = Integer::decrypt(&res[1], &client_key);
    let elasped = now.elapsed();
    println!(
        "{}^-1 % {} -> {}",
        format(x1),
        format(p),
        format(res1_decoded)
    );
    println!("should be {}", format(inverse_mod_native(x1, p)));
    println!(
        "{}^-1 % {} -> {}",
        format(y1),
        format(p),
        format(res2_decoded)
    );
    println!("should be {}", format(inverse_mod_native(y1, p)));
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
    // println!("should be {}", format(pow_mod_native(x1, y1, p)));
    // println!("pow mod in {:.2} s\n", elasped.as_secs_f32());

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

    //let now = Instant::now();
    //let (x_new, y_new, z_new) = group_projective_add_projective::<NUM_BLOCK, _>(
    //&ct_x1,
    //&ct_y1,
    //&client_key.encrypt_radix(1, NUM_BLOCK),
    //&ct_x2,
    //&ct_y2,
    //&client_key.encrypt_radix(1, NUM_BLOCK),
    //p,
    //&server_key,
    //);
    //let elasped = now.elapsed();
    //let x_dec = client_key.decrypt_radix::<Integer>(&x_new);
    //let y_dec = client_key.decrypt_radix::<Integer>(&y_new);
    //let z_dec = client_key.decrypt_radix::<Integer>(&z_new);
    //println!(
    //"{},{},{} + {},{},{} -> {},{},{}",
    //format(x1),
    //format(y1),
    //format(1),
    //format(x2),
    //format(y2),
    //format(1),
    //format(x_dec),
    //format(y_dec),
    //format(z_dec)
    //);
    //println!("group add in {} s", elasped.as_secs_f32());

    //let now = Instant::now();
    //let (x_new, y_new, z_new) = group_projective_scalar_mul_constant_windowed::<8, NUM_BLOCK, _>(
    //x1,
    //y1,
    //&ct_x2,
    //p,
    //&server_key,
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
    //println!("Stats: {}", ProtocolStats::stats());

    //let now = Instant::now();
    //let (x_aff, y_aff) =
    //group_projective_into_affine::<NUM_BLOCK, _>(&x_new, &y_new, &z_new, p, &server_key);
    //let elasped = now.elapsed();
    //println!(
    //"{},{},{} -> {},{}",
    //format(x_dec),
    //format(y_dec),
    //format(z_dec),
    //format(client_key.decrypt_radix::<Integer>(&x_aff)),
    //format(client_key.decrypt_radix::<Integer>(&y_aff))
    //);
    //println!(
    //"group projective into affine in {:.2}s",
    //elasped.as_secs_f32()
    //);
}
