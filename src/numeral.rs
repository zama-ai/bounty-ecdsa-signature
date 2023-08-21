use num_bigint::BigInt;
use tfhe::{
    core_crypto::prelude::{Numeric, UnsignedNumeric},
    integer::{
        block_decomposition::{DecomposableInto, RecomposableFrom},
        server_key::{MiniUnsignedInteger, Reciprocable, ScalarMultiplier, TwosComplementNegation},
        ClientKey, RadixCiphertext, U256,
    },
};

use crate::helper::{format, to_bigint};

pub trait Numeral:
    Numeric
    + UnsignedNumeric
    + DecomposableInto<u64>
    + DecomposableInto<u8>
    + RecomposableFrom<u64>
    + RecomposableFrom<u8>
    + MiniUnsignedInteger
    + Reciprocable
    + ScalarMultiplier
    + TwosComplementNegation
    + Copy
    + Sync
    + Send
{
    fn format(&self) -> String {
        format(*self)
    }

    fn decrypt(ciphertext: &RadixCiphertext, client_key: &ClientKey) -> Self {
        client_key.decrypt_radix::<Self>(ciphertext)
    }

    fn decrypt_bigint(ciphertext: &RadixCiphertext, client_key: &ClientKey) -> BigInt {
        to_bigint(client_key.decrypt_radix::<U256>(ciphertext))
    }
}

impl<T> Numeral for T where
    T: Numeric
        + UnsignedNumeric
        + DecomposableInto<u64>
        + DecomposableInto<u8>
        + RecomposableFrom<u64>
        + RecomposableFrom<u8>
        + MiniUnsignedInteger
        + Reciprocable
        + ScalarMultiplier
        + TwosComplementNegation
        + Copy
        + Sync
        + Send
{
}
