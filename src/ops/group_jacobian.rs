use std::time::Instant;

use rand::Rng;
use tfhe::{
    core_crypto::prelude::Numeric,
    integer::{
        block_decomposition::{DecomposableInto, RecomposableFrom},
        RadixCiphertext, ServerKey,
    },
};

use crate::helper::{format, read_client_key};

use super::{add_mod, double_mod, inverse_mod, mul_mod, square_mod, sub_mod};

/// Projective point at infinity in jacobian coordinates
pub fn group_projective_zero<const NB: usize>(
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext, RadixCiphertext) {
    (
        server_key.create_trivial_radix(1, NB),
        server_key.create_trivial_radix(1, NB),
        server_key.create_trivial_radix(0, NB),
    )
}

pub fn group_projective_double<
    const NB: usize,
    P: DecomposableInto<u64> + RecomposableFrom<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    x: &RadixCiphertext,
    y: &RadixCiphertext,
    z: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext, RadixCiphertext) {
    // case curve a = 0
    // a = x^2
    // b = y^2

    #[cfg(feature = "high_level_timing")]
    let ops_start = Instant::now();
    #[cfg(feature = "high_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);

    #[cfg(feature = "high_level_timing")]
    println!("group projective double jacobian start -- ref {}", task_ref);

    let (a, b) = rayon::join(
        || square_mod::<NB, _>(x, p, server_key),
        || square_mod::<NB, _>(y, p, server_key),
    );
    // c = b^2
    let c = square_mod::<NB, _>(&b, p, server_key);
    // d = 2*((x + b)^2-(a + c))
    let (xb, ac) = rayon::join(
        || add_mod::<NB, _>(x, &b, p, server_key),
        || add_mod::<NB, _>(&a, &c, p, server_key),
    );
    let d = double_mod::<NB, _>(
        &sub_mod::<NB, _>(&square_mod::<NB, _>(&xb, p, server_key), &ac, p, server_key),
        p,
        server_key,
    );
    // e = 3*a
    let e = add_mod::<NB, _>(&double_mod::<NB, _>(&a, p, server_key), &a, p, server_key);
    // f = e^2
    // z' = 2*y*z
    let (f, z_prime) = rayon::join(
        || square_mod::<NB, _>(&e, p, server_key),
        || double_mod::<NB, _>(&mul_mod::<NB, _>(y, z, p, server_key), p, server_key),
    );
    // x' = f - 2*d
    let x_prime = sub_mod::<NB, _>(&f, &double_mod::<NB, _>(&d, p, server_key), p, server_key);
    // y' = e*(d - x') - 8*c
    let (edx, c8) = rayon::join(
        || {
            mul_mod::<NB, _>(
                &e,
                &sub_mod::<NB, _>(&d, &x_prime, p, server_key),
                p,
                server_key,
            )
        },
        || {
            let c2 = double_mod::<NB, _>(&c, p, server_key);
            let c4 = double_mod::<NB, _>(&c2, p, server_key);
            double_mod::<NB, _>(&c4, p, server_key)
        },
    );
    let y_prime = sub_mod::<NB, _>(&edx, &c8, p, server_key);
    #[cfg(feature = "high_level_timing")]
    println!(
        "group projective double done in {:.2}s -- ref {}",
        ops_start.elapsed().as_secs_f64(),
        task_ref
    );
    (x_prime, y_prime, z_prime)
}

#[allow(clippy::too_many_arguments)]
pub fn group_projective_add_projective<
    const NB: usize,
    P: DecomposableInto<u64> + RecomposableFrom<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    x0: &RadixCiphertext,
    y0: &RadixCiphertext,
    z0: &RadixCiphertext,
    x1: &RadixCiphertext,
    y1: &RadixCiphertext,
    z1: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext, RadixCiphertext) {
    #[cfg(feature = "high_level_timing")]
    let ops_start = Instant::now();
    #[cfg(feature = "high_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);
    #[cfg(feature = "high_level_timing")]
    println!("group projective add jacobian start -- ref {}", task_ref);

    // z0z0 = z0^2
    // z1z1 = z1^2
    let (z0z0, z1z1) = rayon::join(
        || square_mod::<NB, _>(z0, p, server_key),
        || square_mod::<NB, _>(z1, p, server_key),
    );
    // u0 = x0*z1z1
    // u1 = x1*z0z0
    let (u0, u1) = rayon::join(
        || mul_mod::<NB, _>(x0, &z1z1, p, server_key),
        || mul_mod::<NB, _>(x1, &z0z0, p, server_key),
    );
    // s0 = y0*z1*z1z1
    // s1 = y1*z0*z0z0
    let (s0, s1) = rayon::join(
        || {
            mul_mod::<NB, _>(
                y0,
                &mul_mod::<NB, _>(z1, &z1z1, p, server_key),
                p,
                server_key,
            )
        },
        || {
            mul_mod::<NB, _>(
                y1,
                &mul_mod::<NB, _>(z0, &z0z0, p, server_key),
                p,
                server_key,
            )
        },
    );
    // h = u1 - u0
    // r = 2*(s1 - s0)
    let (h, r) = rayon::join(
        || sub_mod::<NB, _>(&u1, &u0, p, server_key),
        || double_mod::<NB, _>(&sub_mod::<NB, _>(&s1, &s0, p, server_key), p, server_key),
    );
    // i = (2*h)^2
    let i = square_mod::<NB, _>(&double_mod::<NB, _>(&h, p, server_key), p, server_key);
    // j = h*i
    // v = u0*i
    let (j, v) = rayon::join(
        || mul_mod::<NB, _>(&h, &i, p, server_key),
        || mul_mod::<NB, _>(&u0, &i, p, server_key),
    );
    // x_prime = r^2 - j - 2*v
    // y_prime = r*(v - x_prime) - 2*s0*j
    // z_prime = ((z0 + z1)^2 - z0z0 - z1z1)*h
    let ((r2, s0j2), z0z12) = rayon::join(
        || {
            rayon::join(
                || square_mod::<NB, _>(&r, p, server_key),
                || double_mod::<NB, _>(&mul_mod::<NB, _>(&s0, &j, p, server_key), p, server_key),
            )
        },
        || square_mod::<NB, _>(&add_mod::<NB, _>(z0, z1, p, server_key), p, server_key),
    );
    let mut x_prime = sub_mod::<NB, _>(
        &sub_mod::<NB, _>(&r2, &j, p, server_key),
        &double_mod::<NB, _>(&v, p, server_key),
        p,
        server_key,
    );
    let (mut y_prime, mut z_prime) = rayon::join(
        || {
            sub_mod::<NB, _>(
                &mul_mod::<NB, _>(
                    &r,
                    &sub_mod::<NB, _>(&v, &x_prime, p, server_key),
                    p,
                    server_key,
                ),
                &s0j2,
                p,
                server_key,
            )
        },
        || {
            mul_mod::<NB, _>(
                &sub_mod::<NB, _>(
                    &sub_mod::<NB, _>(&z0z12, &z0z0, p, server_key),
                    &z1z1,
                    p,
                    server_key,
                ),
                &h,
                p,
                server_key,
            )
        },
    );

    ///////////////////////////////// COPY PASTE FROM HOMOGENOUS - might be buggy /////////////////////////////////
    // z1'/z0' 0  1
    //    0    x' x1
    //    1    x0 x0
    // x'' =  x' * is_z0_z1_non_zero + (x0 + x1) * not_is_z0_z1_non_zero
    // y'' =  y' * is_z0_z1_non_zero + (y0 + y1) * not_is_z0_z1_non_zero
    // z'' =  z' * is_z0_z1_non_zero + (z0 + z1) * not_is_z0_z1_non_zero
    let (mut is_z0_non_zero, mut is_z1_non_zero) = rayon::join(
        || server_key.smart_scalar_ne_parallelized(&mut z0.clone(), 0),
        || server_key.smart_scalar_ne_parallelized(&mut z1.clone(), 0),
    );
    server_key.trim_radix_blocks_msb_assign(&mut is_z0_non_zero, NB - 1);
    server_key.trim_radix_blocks_msb_assign(&mut is_z1_non_zero, NB - 1);
    let mut is_z0_z1_non_zero =
        server_key.smart_bitand_parallelized(&mut is_z0_non_zero, &mut is_z1_non_zero);
    let not_is_z0_z1_non_zero = server_key.smart_sub_parallelized(
        &mut server_key.create_trivial_radix(1, 1),
        &mut is_z0_z1_non_zero,
    );

    let (((mut xp1, mut xp2), (mut yp1, mut yp2)), (mut zp1, mut zp2)) = rayon::join(
        || {
            rayon::join(
                || {
                    rayon::join(
                        || {
                            server_key.smart_mul_parallelized(
                                &mut x_prime,
                                &mut is_z0_z1_non_zero.clone(),
                            )
                        },
                        || {
                            server_key.smart_mul_parallelized(
                                &mut server_key
                                    .smart_add_parallelized(&mut x0.clone(), &mut x1.clone()),
                                &mut not_is_z0_z1_non_zero.clone(),
                            )
                        },
                    )
                },
                || {
                    rayon::join(
                        || {
                            server_key.smart_mul_parallelized(
                                &mut y_prime,
                                &mut is_z0_z1_non_zero.clone(),
                            )
                        },
                        || {
                            server_key.smart_mul_parallelized(
                                &mut server_key
                                    .smart_add_parallelized(&mut y0.clone(), &mut y1.clone()),
                                &mut not_is_z0_z1_non_zero.clone(),
                            )
                        },
                    )
                },
            )
        },
        || {
            rayon::join(
                || server_key.smart_mul_parallelized(&mut z_prime, &mut is_z0_z1_non_zero.clone()),
                || {
                    server_key.smart_mul_parallelized(
                        &mut server_key.smart_add_parallelized(&mut z0.clone(), &mut z1.clone()),
                        &mut not_is_z0_z1_non_zero.clone(),
                    )
                },
            )
        },
    );

    ((x_prime, y_prime), z_prime) = rayon::join(
        || {
            rayon::join(
                || server_key.smart_add_parallelized(&mut xp1, &mut xp2),
                || server_key.smart_add_parallelized(&mut yp1, &mut yp2),
            )
        },
        || server_key.smart_add_parallelized(&mut zp1, &mut zp2),
    );
    ///////////////////////////////// END OF COPY PASTE FROM HOMOGENOUS - might be buggy /////////////////////////////////

    #[cfg(feature = "high_level_timing")]
    println!(
        "group projective add done in {:.2}s -- ref {}",
        ops_start.elapsed().as_secs_f64(),
        task_ref
    );

    (x_prime, y_prime, z_prime)
}

pub fn group_projective_scalar_mul<
    const NB: usize,
    P: DecomposableInto<u64> + RecomposableFrom<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    x: &RadixCiphertext,
    y: &RadixCiphertext,
    z: &RadixCiphertext,
    scalar: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext, RadixCiphertext) {
    #[cfg(feature = "high_level_timing")]
    let ops_start = Instant::now();
    #[cfg(feature = "high_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);
    #[cfg(feature = "high_level_timing")]
    println!(
        "group projective scalar mul jacobian start -- ref {}",
        task_ref
    );

    let mut tmp_x = x.clone();
    let mut tmp_y = y.clone();
    let mut tmp_z = z.clone();
    let mut scalar = scalar.clone();
    let mut res_x = server_key.create_trivial_radix(1, NB);
    let mut res_y = server_key.create_trivial_radix(1, NB);
    let mut res_z = server_key.create_trivial_radix(0, NB);

    for _i in 0..<P as Numeric>::BITS {
        #[cfg(feature = "high_level_timing")]
        let bit_start = Instant::now();

        let (mut bit, new_scalar) = rayon::join(
            || server_key.scalar_bitand_parallelized(&scalar, 1),
            || server_key.scalar_right_shift_parallelized(&scalar, 1),
        );
        server_key.trim_radix_blocks_msb_assign(&mut bit, NB - 1);
        scalar = new_scalar;
        ((res_x, res_y, res_z), (tmp_x, tmp_y, tmp_z)) = rayon::join(
            || {
                let ((x_to_add, y_to_add), z_to_add) = rayon::join(
                    || {
                        rayon::join(
                            || {
                                let mut x = server_key
                                    .smart_mul_parallelized(&mut tmp_x.clone(), &mut bit.clone());
                                server_key.scalar_add_assign_parallelized(&mut x, 1);
                                server_key.smart_sub_assign_parallelized(&mut x, &mut bit.clone());
                                x
                            },
                            || {
                                let mut y = server_key
                                    .smart_mul_parallelized(&mut tmp_y.clone(), &mut bit.clone());
                                server_key.scalar_add_assign_parallelized(&mut y, 1);
                                server_key.smart_sub_assign_parallelized(&mut y, &mut bit.clone());
                                y
                            },
                        )
                    },
                    || server_key.smart_mul_parallelized(&mut tmp_z.clone(), &mut bit.clone()),
                );
                group_projective_add_projective::<NB, _>(
                    &res_x, &res_y, &res_z, &x_to_add, &y_to_add, &z_to_add, p, server_key,
                )
            },
            || group_projective_double::<NB, _>(&tmp_x, &tmp_y, &tmp_z, p, server_key),
        );
        #[cfg(feature = "high_level_timing")]
        read_client_key(|client_key| {
            println!("Bit = {}", format(client_key.decrypt_radix::<P>(&bit)),);
            println!(
                "Res {},{},{}",
                format(client_key.decrypt_radix::<P>(&res_x)),
                format(client_key.decrypt_radix::<P>(&res_y)),
                format(client_key.decrypt_radix::<P>(&res_z)),
            );
            println!(
                "Tmp {},{},{}",
                format(client_key.decrypt_radix::<P>(&tmp_x)),
                format(client_key.decrypt_radix::<P>(&tmp_y)),
                format(client_key.decrypt_radix::<P>(&tmp_z)),
            );
            println!(
                "----Scalar mul bit {_i} done in {:.2}s -- ref {}",
                bit_start.elapsed().as_secs_f32(),
                task_ref
            );
        });
    }

    #[cfg(feature = "high_level_timing")]
    println!(
        "group projective scalar mul done in {:.2}s -- ref {}",
        ops_start.elapsed().as_secs_f64(),
        task_ref
    );
    (res_x, res_y, res_z)
}

//pub fn group_projective_scalar_mul_constant<
//const NB: usize,
//P: DecomposableInto<u64> + RecomposableFrom<u64> + DecomposableInto<u8> + Copy + Sync,
//>(
//x: P,
//y: P,
//z: P,
//scalar: &RadixCiphertext,
//p: P,
//server_key: &ServerKey,
//) -> (RadixCiphertext, RadixCiphertext, RadixCiphertext) {
//#[cfg(feature = "high_level_timing")]
//let ops_start = Instant::now();
//#[cfg(feature = "high_level_timing")]
//let task_ref = rand::thread_rng().gen_range(0..1000);
//#[cfg(feature = "high_level_timing")]
//println!(
//"group projective scalar mul jacobian start -- ref {}",
//task_ref
//);

//let mut tmp_x = x;
//let mut tmp_y = y;
//let mut tmp_z = z;
//let mut scalar = scalar.clone();
//let mut res_x = server_key.create_trivial_radix(1, NB);
//let mut res_y = server_key.create_trivial_radix(1, NB);
//let mut res_z = server_key.create_trivial_radix(0, NB);

//for _i in 0..<P as Numeric>::BITS {
//#[cfg(feature = "high_level_timing")]
//let bit_start = Instant::now();

//let (mut bit, new_scalar) = rayon::join(
//|| server_key.scalar_bitand_parallelized(&scalar, 1),
//|| server_key.scalar_right_shift_parallelized(&scalar, 1),
//);
//server_key.trim_radix_blocks_msb_assign(&mut bit, NB - 1);
//scalar = new_scalar;
//let ((x_to_add, y_to_add), z_to_add) = rayon::join(
//|| {
//rayon::join(
//|| {
//let mut x =
//server_key.smart_mul_parallelized(&mut tmp_x.clone(), &mut bit.clone());
//server_key.scalar_add_assign_parallelized(&mut x, 1);
//server_key.smart_sub_assign_parallelized(&mut x, &mut bit.clone());
//x
//},
//|| {
//let mut y =
//server_key.smart_mul_parallelized(&mut tmp_y.clone(), &mut bit.clone());
//server_key.scalar_add_assign_parallelized(&mut y, 1);
//server_key.smart_sub_assign_parallelized(&mut y, &mut bit.clone());
//y
//},
//)
//},
//|| server_key.smart_mul_parallelized(&mut tmp_z.clone(), &mut bit.clone()),
//);
//(res_x, res_y, res_z) = group_projective_add_projective::<NB, _>(
//&res_x, &res_y, &res_z, &x_to_add, &y_to_add, &z_to_add, p, server_key,
//);
//#[cfg(feature = "high_level_timing")]
//read_client_key(|client_key| {
//println!("Bit = {}", format(client_key.decrypt_radix::<P>(&bit)),);
//println!(
//"Res {},{},{}",
//format(client_key.decrypt_radix::<P>(&res_x)),
//format(client_key.decrypt_radix::<P>(&res_y)),
//format(client_key.decrypt_radix::<P>(&res_z)),
//);
//println!("Tmp {},{},{}", format(tmp_x), format(tmp_y), format(tmp_z),);
//println!(
//"----Scalar mul bit {_i} done in {:.2}s -- ref {}",
//bit_start.elapsed().as_secs_f32(),
//task_ref
//);
//});
//}

//#[cfg(feature = "high_level_timing")]
//println!(
//"group projective scalar mul done in {:.2}s -- ref {}",
//ops_start.elapsed().as_secs_f64(),
//task_ref
//);
//(res_x, res_y, res_z)
//}

pub fn group_projective_into_affine<
    const NB: usize,
    P: DecomposableInto<u64> + RecomposableFrom<u64> + DecomposableInto<u8> + Copy + Sync,
>(
    x: &RadixCiphertext,
    y: &RadixCiphertext,
    z: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> (RadixCiphertext, RadixCiphertext) {
    #[cfg(feature = "high_level_timing")]
    let ops_start = Instant::now();
    #[cfg(feature = "high_level_timing")]
    let task_ref = rand::thread_rng().gen_range(0..1000);
    #[cfg(feature = "high_level_timing")]
    println!(
        "group projective into affine jacobian start -- ref {}",
        task_ref
    );

    let z_inv = inverse_mod::<NB, _>(z, p, server_key);
    let z_inv2 = square_mod::<NB, _>(&z_inv, p, server_key);
    let z_inv3 = mul_mod::<NB, _>(&z_inv2, &z_inv, p, server_key);

    let res = rayon::join(
        || mul_mod::<NB, _>(x, &z_inv2, p, server_key),
        || mul_mod::<NB, _>(y, &z_inv3, p, server_key),
    );
    #[cfg(feature = "high_level_timing")]
    println!(
        "group projective into affine done in {:.2}s -- ref {}",
        ops_start.elapsed().as_secs_f64(),
        task_ref
    );
    res
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use tfhe::{integer::keycache::IntegerKeyCache, shortint::prelude::PARAM_MESSAGE_2_CARRY_2};

    use crate::helper::{format, set_client_key};

    use super::{
        group_projective_add_projective, group_projective_double, group_projective_into_affine,
        group_projective_scalar_mul,
    };

    #[test]
    fn correct_jacobian_double() {
        let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);

        const NUM_BLOCK: usize = 4;
        type Integer = u8;
        let p: Integer = 251;
        let x1: Integer = 8;
        let y1: Integer = 45;
        //let x2: Integer = 26;
        //let y2: Integer = 55;

        let ct_x1 = client_key.encrypt_radix(x1, NUM_BLOCK);
        let ct_y1 = client_key.encrypt_radix(y1, NUM_BLOCK);
        //let ct_x2 = client_key.encrypt_radix(x2, NUM_BLOCK);
        //let ct_y2 = client_key.encrypt_radix(y2, NUM_BLOCK);

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

    #[test]
    fn correct_jacobian_add() {
        let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);

        const NUM_BLOCK: usize = 4;
        type Integer = u8;
        let p: Integer = 251;
        let x1: Integer = 8;
        let y1: Integer = 45;
        let x2: Integer = 26;
        let y2: Integer = 55;

        let ct_x1 = client_key.encrypt_radix(x1, NUM_BLOCK);
        let ct_y1 = client_key.encrypt_radix(y1, NUM_BLOCK);
        let ct_x2 = client_key.encrypt_radix(x2, NUM_BLOCK);
        let ct_y2 = client_key.encrypt_radix(y2, NUM_BLOCK);

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

    #[test]
    fn correct_jacobian_scalar_mul() {
        let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);
        set_client_key(&client_key);

        const NUM_BLOCK: usize = 4;
        type Integer = u8;
        let p: Integer = 251;
        let x1: Integer = 8;
        let y1: Integer = 45;
        let x2: Integer = 26;

        let ct_x1 = client_key.encrypt_radix(x1, NUM_BLOCK);
        let ct_y1 = client_key.encrypt_radix(y1, NUM_BLOCK);
        let ct_x2 = client_key.encrypt_radix(x2, NUM_BLOCK);

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
}
