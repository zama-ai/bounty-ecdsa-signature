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
    use crate::ops::native::{add_mod_native, inverse_mod_native, sub_mod_native};

    use super::mul_mod_native;

    #[test]
    fn correct_native_op() {
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
    }
}
