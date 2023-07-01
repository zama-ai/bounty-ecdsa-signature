use num_bigint::BigInt;
use tfhe::integer::{
    block_decomposition::{BlockDecomposer, DecomposableInto},
    U256,
};

pub fn format<T: DecomposableInto<u8> + Copy>(a: T) -> String {
    BigInt::from_bytes_le(
        num_bigint::Sign::Plus,
        &BlockDecomposer::new(a, 8)
            .iter_as::<u8>()
            .collect::<Vec<_>>(),
    )
    .to_string()
}
