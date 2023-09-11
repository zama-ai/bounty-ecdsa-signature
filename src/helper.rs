use num_bigint::BigInt;
use tfhe::integer::{
    block_decomposition::{BlockDecomposer, DecomposableInto, RecomposableFrom},
    ClientKey, U256,
};

use crate::CLIENT_KEY;

pub fn bigint_ilog2_ceil(value: &BigInt) -> u32 {
    let mut value = value.clone();
    let mut i = 0;
    if value == BigInt::from(0) {
        0
    } else if value == BigInt::from(1) {
        1
    } else {
        while value > BigInt::from(0) {
            value >>= 1;
            i += 1;
        }
        i
    }
}

pub fn bigint_ilog2_floor(value: &BigInt) -> u32 {
    let mut value = value.clone();
    let mut i = 0;
    if value == BigInt::from(0) {
        0
    } else if value == BigInt::from(1) {
        1
    } else {
        while value > BigInt::from(1) {
            value >>= 1;
            i += 1;
        }
        i
    }
}

pub fn format<T: DecomposableInto<u8> + Copy>(a: T) -> String {
    BigInt::from_bytes_le(
        num_bigint::Sign::Plus,
        &BlockDecomposer::new(a, 8)
            .iter_as::<u8>()
            .collect::<Vec<_>>(),
    )
    .to_string()
}

pub fn to_bigint<T: DecomposableInto<u8>>(a: T) -> BigInt {
    BigInt::from_bytes_le(
        num_bigint::Sign::Plus,
        &BlockDecomposer::new(a, 8)
            .iter_as::<u8>()
            .collect::<Vec<_>>(),
    )
}

pub fn from_bigint<T: DecomposableInto<u8> + RecomposableFrom<u8>>(a: &BigInt) -> T {
    let mut res = T::ZERO;
    for (i, b) in a.to_bytes_le().1.iter().enumerate() {
        res += T::cast_from(*b) << (i * 8) as u32;
    }
    res
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
        res += U256::from(*b) << (i * 8) as u32;
    }
    res
}

pub fn bigint_to_u128(a: &BigInt) -> u128 {
    let mut res = 0;
    for (i, b) in a.to_bytes_le().1.iter().enumerate() {
        res += (*b as u128) << (i * 8);
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
    if let Some(client_key) = CLIENT_KEY.read().unwrap().as_ref() {
        f(client_key);
    }
}

pub fn set_client_key(client_key: &ClientKey) {
    *CLIENT_KEY.write().unwrap() = Some(client_key.clone());
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use num_bigint::BigInt;

    use crate::helper::{bigint_ilog2_ceil, bigint_ilog2_floor, format};

    use super::u256_from_decimal_string;

    #[test]
    fn correct_u256_from_string() {
        let str_value = "47328345983454384985382623486293476776767";
        let u256_value = u256_from_decimal_string(str_value);
        assert_eq!(format(u256_value), str_value);
    }

    #[test]
    fn correct_bigint_ilog2() {
        let value = 1237582375u128;
        let bigint_value = BigInt::from(value);
        let res = bigint_ilog2_floor(&bigint_value);
        assert_eq!(res, value.ilog2());
        assert_eq!(bigint_ilog2_ceil(&BigInt::from_str("115792089237316195423570985008687907853269984665640564039457584007908834671663").unwrap()), 256);
    }
}
