use std::time::Instant;

use tfhe::integer::ciphertext::BaseRadixCiphertext;
use tfhe::integer::keycache::IntegerKeyCache;
use tfhe::integer::{gen_keys_radix, ServerKey};
use tfhe::shortint::parameters::PARAM_MESSAGE_2_CARRY_2;
use tfhe::shortint::Ciphertext;

type Ctxt = BaseRadixCiphertext<Ciphertext>;

fn add_mod(a: &Ctxt, b: &Ctxt, p: &Ctxt, server_key: &ServerKey) -> Ctxt {
    // assume large p and a,b < p
    let added = server_key.add_parallelized(&a, &b);
    let is_gt = server_key.gt_parallelized(&added, &p);
    let to_sub = server_key.mul_parallelized(&p, &is_gt);
    let res = server_key.sub_parallelized(&added, &to_sub);
    res
}

fn mul_mod(a: &Ctxt, b: &Ctxt, p: &Ctxt, server_key: &ServerKey) -> Ctxt {
    // assume large p and a,b < p
    let mut res = server_key.create_trivial_radix(0u64, 4);
    let mut tmp = a.clone();
    let mut bb = b.clone();
    let now = Instant::now();

    for i in 0..8 {
        let bit = server_key.scalar_bitand_parallelized(&bb, 1);
        let to_add = server_key.mul_parallelized(&tmp, &bit);

        let mut result_res = None;
        let mut result_tmp = None;
        let mut result_bb = None;

        rayon::scope(|s| {
            s.spawn(|_| result_res = Some(add_mod(&res, &to_add, &p, &server_key)));
            s.spawn(|_| result_tmp = Some(add_mod(&tmp, &tmp, &p, &server_key)));
            s.spawn(|_| result_bb = Some(server_key.scalar_right_shift_parallelized(&bb, 1)));
        });

        res = result_res.unwrap();
        tmp = result_tmp.unwrap();
        bb = result_bb.unwrap();

        println!("time used {i} - {}", now.elapsed().as_secs());
    }
    res
}

fn main() {
    let num_block = 4;
    let (client_key, server_key) = IntegerKeyCache.get_from_params(PARAM_MESSAGE_2_CARRY_2);
    // let (client_key, server_key) = gen_keys_radix(&PARAM_MESSAGE_2_CARRY_2, num_block);
    let msg1 = 30u64;
    let msg2 = 159u64;
    let p = 191u64;

    // We use the client key to encrypt two messages:
    let ct_1 = client_key.encrypt_radix(msg1, num_block);
    let ct_2 = client_key.encrypt_radix(msg2, num_block);
    let ct_p = server_key.create_trivial_radix(p, num_block);

    // let now = Instant::now();
    // let res = add_mod(&ct_1, &ct_2, &ct_p, &server_key);
    // let elasped = now.elapsed();
    // let ddd: u64 = client_key.decrypt(&res);
    // print!("{} + {} mod {} -> {}", msg1, msg2, p, ddd);
    // print!("in {} ms\n", elasped.as_secs());

    let now = Instant::now();
    let res = mul_mod(&ct_1, &ct_2, &ct_p, &server_key);
    let elasped = now.elapsed();
    let ddd: u64 = client_key.decrypt_radix(&res);
    print!("{} * {} mod {} -> {}", msg1, msg2, p, ddd);
    print!("in {} ms\n", elasped.as_secs());
}
