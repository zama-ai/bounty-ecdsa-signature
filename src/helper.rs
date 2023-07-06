use num_bigint::BigInt;
use tfhe::integer::{
    block_decomposition::{BlockDecomposer, DecomposableInto},
    ClientKey, U256,
};

use crate::CLIENT_KEY;

pub fn format<T: DecomposableInto<u8> + Copy>(a: T) -> String {
    BigInt::from_bytes_le(
        num_bigint::Sign::Plus,
        &BlockDecomposer::new(a, 8)
            .iter_as::<u8>()
            .collect::<Vec<_>>(),
    )
    .to_string()
}

pub fn u256_to_bigint(a: U256) -> BigInt {
    BigInt::from_bytes_le(
        num_bigint::Sign::Plus,
        &BlockDecomposer::new(a, 8)
            .iter_as::<u8>()
            .collect::<Vec<_>>(),
    )
}

pub fn bigint_to_u256(a: &BigInt) -> U256 {
    let mut res = U256::ZERO;
    for (i, b) in a.to_bytes_le().1.iter().enumerate() {
        res += U256::from(*b as u8) << (i * 8) as u32;
    }
    res
}

pub fn u256_from_decimal_string(s: &str) -> U256 {
    let mut res = U256::ZERO;
    for c in s.chars() {
        res *= 10u8.into();
        res += (c as u8 - b'0').into();
    }
    res
}

pub fn read_client_key<F: FnOnce(&ClientKey)>(f: F) {
    CLIENT_KEY
        .read()
        .unwrap()
        .as_ref()
        .inspect(|client_key| f(client_key));
}

#[cfg(test)]
mod tests {
    use crate::helper::format;

    use super::u256_from_decimal_string;

    #[test]
    fn correct_u256_from_string() {
        let str_value = "47328345983454384985382623486293476776767";
        let u256_value = u256_from_decimal_string(str_value);
        assert_eq!(format(u256_value), str_value);
    }
}
