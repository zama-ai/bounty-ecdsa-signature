use num_bigint::BigInt;
use tfhe::integer::block_decomposition::{DecomposableInto, RecomposableFrom};

use crate::{
    helper::{from_bigint, to_bigint},
    numeral::Numeral,
};

/// a^-1 mod p where a*a^-1 = 1 mod p
#[inline]
pub fn inverse_mod_native<P: Numeral>(a: P, p: P) -> P {
    pow_mod_native(a, p - P::TWO, p)
}

#[inline]
pub fn modulo_native<P: Numeral>(a: P, p: P) -> P {
    let a_bigint = to_bigint(a);
    let p_bigint = to_bigint(p);
    from_bigint(&(&a_bigint % &p_bigint))
}

/// a^b mod p
#[inline]
pub fn pow_mod_native<P: Numeral>(a: P, b: P, p: P) -> P {
    let a_bigint = to_bigint(a);
    let b_bigint = to_bigint(b);
    let p_bigint = to_bigint(p);

    let mut res = BigInt::from(1);
    let mut base = &a_bigint % &p_bigint;
    let mut exponent = b_bigint;

    while exponent > BigInt::from(0) {
        if &exponent % 2 == BigInt::from(1) {
            res = (&res * &base) % &p_bigint;
        }
        exponent >>= 1;
        base = (&base * &base) % &p_bigint;
    }

    from_bigint(&res)
}

/// a + b mod p
pub fn add_mod_native<P: Numeral>(a: P, b: P, p: P) -> P {
    let a_bigint = to_bigint(a);
    let b_bigint = to_bigint(b);
    let p_bigint = to_bigint(p);
    from_bigint(&((a_bigint + b_bigint) % p_bigint))
}

/// a - b mod p
pub fn sub_mod_native<P: Numeral>(a: P, b: P, p: P) -> P {
    let a_bigint = to_bigint(a);
    let b_bigint = to_bigint(b);
    let p_bigint = to_bigint(p);
    if a_bigint < b_bigint {
        from_bigint(&((a_bigint + p_bigint) - b_bigint))
    } else {
        from_bigint(&((a_bigint - b_bigint) % p_bigint))
    }
}

/// a * b mod p
pub fn mul_mod_native<P: Numeral>(a: P, b: P, p: P) -> P {
    let a_bigint = to_bigint(a);
    let b_bigint = to_bigint(b);
    let p_bigint = to_bigint(p);
    from_bigint(&((a_bigint * b_bigint) % p_bigint))
}

/// a^2 mod p
pub fn square_mod_native<P: Numeral>(a: P, p: P) -> P {
    let a_bigint = to_bigint(a);
    let p_bigint = to_bigint(p);
    from_bigint(&((&a_bigint * &a_bigint) % p_bigint))
}

/// a*2 mod p
pub fn double_mod_native<P: Numeral>(a: P, p: P) -> P {
    let a_bigint = to_bigint(a);
    let p_bigint = to_bigint(p);
    from_bigint(&((&a_bigint * 2) % p_bigint))
}

#[cfg(test)]
mod tests {
    use tfhe::integer::U256;

    use crate::{
        helper::u256_from_decimal_string,
        ops::{
            native::{
                add_mod_native, double_mod_native, inverse_mod_native, pow_mod_native,
                square_mod_native, sub_mod_native,
            },
            secp256k1::prelude::FQ_MODULO,
        },
    };

    use super::mul_mod_native;

    #[test]
    fn correct_native_op_small() {
        let a: u8 = 123;
        let b: u8 = 234;
        let p: u8 = 251;

        let res = mul_mod_native(a, b, p);
        assert_eq!(res, 168);

        let res = add_mod_native(a, b, p);
        assert_eq!(res, 106);

        let res = sub_mod_native(a, b, p);
        assert_eq!(res, 140);

        let res = inverse_mod_native(a, p);
        assert_eq!(mul_mod_native(res, a, p), 1);

        let res = double_mod_native(a, p);
        assert_eq!(res, 246);

        let res = square_mod_native(a, p);
        assert_eq!(res, 69);

        let res = pow_mod_native(a, b, p);
        assert_eq!(res, 91);
    }

    #[test]
    fn correct_native_op_big() {
        let p: U256 = *FQ_MODULO;
        let a: U256 = u256_from_decimal_string(
            "158972629851468960855479098042189567798917817837573660423710583832714848",
        );
        let b: U256 = u256_from_decimal_string(
            "65108744961846543415519418389643270459525907322081164366671650776835723265410",
        );

        let res = mul_mod_native(a, b, p);
        assert_eq!(
            res,
            u256_from_decimal_string(
                "114526699765647557269583289524217720249718105067902256167975611404180019746063"
            )
        );

        let res = add_mod_native(a, b, p);
        assert_eq!(
            res,
            u256_from_decimal_string(
                "65108903934476394884480273868741312649093706239899001940332074487419555980258"
            )
        );

        let res = sub_mod_native(a, b, p);
        assert_eq!(
            res,
            u256_from_decimal_string(
                "50683503248099503477012422098142679583311876261377237246446356941656944121101"
            )
        );

        let res = inverse_mod_native(a, p);
        assert_eq!(mul_mod_native(res, a, p), u256_from_decimal_string("1"));

        let res = double_mod_native(a, p);
        assert_eq!(
            res,
            u256_from_decimal_string(
                "317945259702937921710958196084379135597835635675147320847421167665429696"
            )
        );

        let res = square_mod_native(a, p);
        assert_eq!(
            res,
            u256_from_decimal_string(
                "91281256522923538041967617230875549618244827250923126868698688129566942900200"
            )
        );

        let res = pow_mod_native(a, b, p);
        assert_eq!(
            res,
            u256_from_decimal_string(
                "10676208167660690292400660450238079909328736345970444101577880711570166245327"
            )
        );
    }
}
