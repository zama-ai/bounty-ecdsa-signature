pub mod prelude {
    use lazy_static::lazy_static;
    use tfhe::integer::U256;

    use crate::helper::u256_from_decimal_string;

    lazy_static! {
        /// The base prime field modulus of secp256k1
        pub static ref FQ_MODULO: U256 = u256_from_decimal_string("115792089237316195423570985008687907853269984665640564039457584007908834671663");
        /// The scalar prime field modulus of secp256k1
        pub static ref FR_MODULO: U256 = u256_from_decimal_string("115792089237316195423570985008687907852837564279074904382605163141518161494337");
        /// The prime subgroup generator of secp256k1
        pub static ref GENERATOR: (U256,U256) = (
            u256_from_decimal_string("55066263022277343669578718895168534326250603453777594175500187360389116729240"),
            u256_from_decimal_string("32670510020758816978083085130507043184471273380659243275938904335757337482424"),
        );
    }
}
