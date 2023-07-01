use std::time::Instant;

use tfhe::{
    core_crypto::prelude::Numeric,
    integer::{block_decomposition::DecomposableInto, RadixCiphertext, ServerKey},
};

pub fn add_mod<const NB: usize, P: DecomposableInto<u64> + DecomposableInto<u8> + Copy + Sync>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    // assume large p and a,b < p
    let mut added = server_key.smart_add_parallelized(&mut a.clone(), &mut b.clone());
    let mut is_gt = server_key.smart_scalar_gt_parallelized(&mut added, p);
    let mut to_sub =
        server_key.smart_mul_parallelized(&mut is_gt, &mut server_key.create_trivial_radix(p, NB));
    let res = server_key.smart_sub_parallelized(&mut added, &mut to_sub);
    res
}

pub fn mul_mod<const NB: usize, P: DecomposableInto<u64> + DecomposableInto<u8> + Copy + Sync>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    // assume large p and a,b < p
    let mut res = server_key.create_trivial_radix(0u64, NB);
    let mut tmp = a.clone();
    let mut bb = b.clone();
    let now = Instant::now();

    for i in 0..<P as Numeric>::BITS {
        ((bb, tmp), res) = rayon::join(
            || {
                rayon::join(
                    || server_key.scalar_right_shift_parallelized(&bb, 1),
                    || add_mod::<NB, P>(&tmp, &tmp, p, &server_key),
                )
            },
            || {
                let mut bit = server_key.smart_scalar_bitand_parallelized(&mut bb.clone(), 1);
                let to_add = server_key.smart_mul_parallelized(&mut tmp.clone(), &mut bit);
                add_mod::<NB, P>(&res, &to_add, p, &server_key)
            },
        );

        println!("time used {i} - {}", now.elapsed().as_secs());
    }
    res
}
