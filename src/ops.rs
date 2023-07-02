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
    server_key.trim_radix_blocks_msb_assign(&mut is_gt, NB - 1);
    let mut to_sub =
        server_key.smart_mul_parallelized(&mut server_key.create_trivial_radix(p, NB), &mut is_gt);

    server_key.smart_sub_parallelized(&mut added, &mut to_sub)
}

pub fn mul_mod<const NB: usize, P: DecomposableInto<u64> + DecomposableInto<u8> + Copy + Sync>(
    a: &RadixCiphertext,
    b: &RadixCiphertext,
    p: P,
    server_key: &ServerKey,
) -> RadixCiphertext {
    // assume large p and a,b < p
    let now = Instant::now();

    let mut res = server_key.create_trivial_radix(0u64, NB);
    let mut a_tmp = a.clone();
    let b_tmp = b.clone();
    let (mut b_next_tmp, mut bit) = rayon::join(
        || server_key.scalar_right_shift_parallelized(&b_tmp, 1),
        || server_key.smart_scalar_bitand_parallelized(&mut b_tmp.clone(), 1),
    );
    let mut to_add_later = res.clone();

    println!("initial cost - {} ms", now.elapsed().as_millis());

    for i in 0..<P as Numeric>::BITS {
        let now = Instant::now();

        ((b_next_tmp, a_tmp), (bit, (res, to_add_later))) = rayon::join(
            || {
                rayon::join(
                    || server_key.scalar_right_shift_parallelized(&b_next_tmp, 1),
                    || add_mod::<NB, _>(&a_tmp, &a_tmp, p, server_key),
                )
            },
            || {
                rayon::join(
                    || server_key.smart_scalar_bitand_parallelized(&mut b_next_tmp.clone(), 1),
                    || {
                        rayon::join(
                            || add_mod::<NB, _>(&res, &to_add_later, p, server_key),
                            || server_key.smart_mul_parallelized(&mut a_tmp.clone(), &mut bit),
                        )
                    },
                )
            },
        );

        println!("time used {i} - {}", now.elapsed().as_secs());
    }

    add_mod::<NB, _>(&res, &to_add_later, p, server_key)
}
