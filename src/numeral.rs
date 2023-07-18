use tfhe::{
    core_crypto::prelude::Numeric,
    integer::{
        block_decomposition::{DecomposableInto, RecomposableFrom},
        ClientKey, RadixCiphertext,
    },
};

use crate::helper::format;

pub trait Numeral:
    Numeric
    + DecomposableInto<u64>
    + DecomposableInto<u8>
    + RecomposableFrom<u64>
    + RecomposableFrom<u8>
    //+ ScalarMultiplier
    //+ TwosComplementNegation
    + Copy
    + Sync
    + Send
{
    fn format(&self) -> String {
        format(*self)
    }

    fn decrypt(ciphertext: &RadixCiphertext, client_key: &ClientKey) -> Self {
        client_key.decrypt_radix::<Self>(&ciphertext)
    }
}

impl<T> Numeral for T where
    T: Numeric
        + DecomposableInto<u64>
        + DecomposableInto<u8>
        + RecomposableFrom<u64>
        + RecomposableFrom<u8>
        //+ ScalarMultiplier
        //+ TwosComplementNegation
        + Copy
        + Sync
        + Send
{
}
