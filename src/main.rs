use std::time::Instant;

use rand::{thread_rng, Rng};
use tfhe::{
    integer::{ciphertext::BaseRadixCiphertext, keycache::IntegerKeyCache, ServerKey, U256},
    shortint::{parameters::PARAM_MESSAGE_2_CARRY_2, Ciphertext},
};

type Ctxt = BaseRadixCiphertext<Ciphertext>;

const NUM_BLOCK: usize = 129;
const BITS: usize = 256;

fn format_u256(a: &U256) -> String {
    let mut bytes = [0u8; 32];
    a.copy_to_be_byte_slice(&mut bytes);
    bytes
        .iter()
        .map(|b| format!("{:02}", b))
        .collect::<Vec<_>>()
        .join("")
}

fn add_mod(a: &Ctxt, b: &Ctxt, p: U256, server_key: &ServerKey) -> Ctxt {
    // assume large p and a,b < p
    let mut added = server_key.smart_add_parallelized(&mut a.clone(), &mut b.clone());
    let mut is_gt = server_key.smart_scalar_gt_parallelized(&mut added, p);
    let mut to_sub = server_key.smart_mul_parallelized(
        &mut is_gt,
        &mut server_key.create_trivial_radix(p, NUM_BLOCK),
    );
    let res = server_key.smart_sub_parallelized(&mut added, &mut to_sub);
    res
}
fn mul_mod(a: &Ctxt, b: &Ctxt, p: U256, server_key: &ServerKey) -> Ctxt {
    // assume large p and a,b < p
    let mut res = server_key.create_trivial_radix(0u64, NUM_BLOCK);
    let mut tmp = a.clone();
    let mut bb = b.clone();
    let now = Instant::now();

    for i in 0..BITS {
        //let mut result_res = None;
        //let mut result_tmp = None;
        //let mut result_bb = None;

        ((bb, tmp), res) = rayon::join(
            || {
                rayon::join(
                    || server_key.scalar_right_shift_parallelized(&bb, 1),
                    || add_mod(&tmp, &tmp, p, &server_key),
                )
            },
            || {
                let mut bit = server_key.smart_scalar_bitand_parallelized(&mut bb.clone(), 1);
                let to_add = server_key.smart_mul_parallelized(&mut tmp.clone(), &mut bit);
                add_mod(&res, &to_add, p, &server_key)
            },
        );
        //rayon::scope(|s| {
        //s.spawn(|_| result_res = Some(add_mod(&res, &to_add, &p, &server_key)));
        //s.spawn(|_| result_tmp = Some(add_mod(&tmp, &tmp, &p, &server_key)));
        //s.spawn(|_| result_bb = Some(server_key.scalar_right_shift_parallelized(&bb, 1)));
        //});

        //res = result_res.unwrap();
        //tmp = result_tmp.unwrap();
        //bb = result_bb.unwrap();

        println!("time used {i} - {}", now.elapsed().as_secs());
    }
    res
}

fn mul_mod_pipe(a: &Ctxt, b: &Ctxt, p: U256, server_key: &ServerKey) -> Ctxt {
    // assume large p and a,b < p
    let now = Instant::now();

    let mut res = server_key.create_trivial_radix(0u64, NUM_BLOCK);
    let mut a_tmp = a.clone();
    let b_tmp = b.clone();
    let mut b_next_tmp = server_key.scalar_right_shift_parallelized(&b_tmp, 1);
    let mut bit = server_key.smart_scalar_bitand_parallelized(&mut b_tmp.clone(), 1);
    let mut bit_next: BaseRadixCiphertext<Ciphertext>;
    let mut to_add_later = server_key.create_trivial_radix(0u64, NUM_BLOCK);
    println!("initial cost - {}ms", now.elapsed().as_millis());

    for i in 0..BITS {
        let now = Instant::now();

        ((b_next_tmp, a_tmp), (bit_next, (res, to_add_later))) = rayon::join(
            || {
                rayon::join(
                    || server_key.scalar_right_shift_parallelized(&b_next_tmp, 1),
                    || add_mod(&a_tmp, &a_tmp, p, &server_key),
                )
            },
            || {
                rayon::join(
                    || server_key.smart_scalar_bitand_parallelized(&mut b_next_tmp.clone(), 1),
                    || {
                        rayon::join(
                            || add_mod(&res, &to_add_later, p, &server_key),
                            || server_key.smart_mul_parallelized(&mut a_tmp.clone(), &mut bit),
                        )
                    },
                )
            },
        );
        bit = bit_next.clone();
        println!("time used {i} - {}", now.elapsed().as_secs());
    }
    res = add_mod(&res, &to_add_later, p, &server_key);
    res
}

fn main() {
    let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);
    // let (client_key, server_key) = gen_keys_radix(&PARAM_MESSAGE_2_CARRY_2, num_block);
    let mut msg1 = U256::ZERO;
    msg1.copy_from_le_byte_slice(&thread_rng().gen::<[u8; BITS / 8]>());
    //let msg2 = 159u64;
    let mut p = U256::ZERO;
    p.copy_from_be_byte_slice(
        &hex::decode("fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f").unwrap(),
    );

    // We use the client key to encrypt two messages:
    let ct_1 = client_key.encrypt_radix(msg1, NUM_BLOCK);
    //let mut ct_2 = client_key.encrypt_radix(msg2, num_block);
    //let ct_p = server_key.create_trivial_radix(p, NUM_BLOCK);

    //let now = Instant::now();
    //let res = add_mod(&ct_1, &ct_1, p, &server_key);
    //let elasped = now.elapsed();
    //let res = client_key.decrypt_radix::<U256>(&res);
    //println!(
    //"{} + {} mod {} -> {}",
    //format_u256(&msg1),
    //format_u256(&msg1),
    //format_u256(&p),
    //format_u256(&res)
    //);
    //println!("in {} s\n", elasped.as_secs());

    let now = Instant::now();
    let res = mul_mod_pipe(&ct_1, &ct_1, p, &server_key);

    let elasped = now.elapsed();
    let res = client_key.decrypt_radix::<U256>(&res);
    print!(
        "{} * {} mod (new) {} -> {}",
        format_u256(&msg1),
        format_u256(&msg1),
        format_u256(&p),
        format_u256(&res)
    );
    println!("in {} s\n", elasped.as_secs());
    let now = Instant::now();
    let res = mul_mod(&ct_1, &ct_1, p, &server_key);

    let elasped = now.elapsed();
    let res = client_key.decrypt_radix::<U256>(&res);
    print!(
        "{} * {} mod {} -> {}",
        format_u256(&msg1),
        format_u256(&msg1),
        format_u256(&p),
        format_u256(&res)
    );
    println!("in {} s\n", elasped.as_secs());
}
