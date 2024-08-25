#![allow(non_snake_case)]

/// Ref: 4.3.1 Parameters of [On Proving Pairings](https://eprint.iacr.org/2024/640.pdf)
use crate::params;
use ark_bn254::{Fq, Fq2};
use ark_ff::Field;
use num_bigint::BigUint;
use num_traits::{FromPrimitive, Num, One, Pow};
use once_cell::sync::Lazy;
use std::clone::Clone;
use std::ops::{Add, Deref, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};
use std::str::FromStr;

// constant modulus of Fq
pub const MODULUS_STR: &str = "30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47";
pub const MODULUS: Lazy<BigUint> = Lazy::new(|| BigUint::from_str_radix(MODULUS_STR, 16).unwrap());

// const X: &'static [u64] = &[4965661367192848881]. See more on: Config::X
pub static X: Lazy<BigUint> = Lazy::new(|| BigUint::from_i128(4965661367192848881).unwrap());
pub static R: Lazy<BigUint> = Lazy::new(|| {
    BigUint::from_str(
        "21888242871839275222246405745257275088548364400416034343698204186575808495617",
    )
    .unwrap()
});
pub static R_INV: Lazy<BigUint> = Lazy::new(|| {
    BigUint::from_str(
        "495819184011867778744231927046742333492451180917315223017345540833046880485481720031136878341141903241966521818658471092566752321606779256340158678675679238405722886654128392203338228575623261160538734808887996935946888297414610216445334190959815200956855428635568184508263913274453942864817234480763055154719338281461936129150171789463489422401982681230261920147923652438266934726901346095892093443898852488218812468761027620988447655860644584419583586883569984588067403598284748297179498734419889699245081714359110559679136004228878808158639412436468707589339209058958785568729925402190575720856279605832146553573981587948304340677613460685405477047119496887534881410757668344088436651291444274840864486870663164657544390995506448087189408281061890434467956047582679858345583941396130713046072603335601764495918026585155498301896749919393",
    )
        .unwrap()
});

pub const H: Lazy<BigUint> = Lazy::new(||
// h = (p^12 - 1) / r
MODULUS.deref().pow(12).sub(BigUint::one()).div(R.clone()));

pub const D: Lazy<BigUint> = Lazy::new(||
// d = gcd(m, h) = 3
BigUint::from_i8(3).unwrap());
pub const M: Lazy<BigUint> = Lazy::new(||
// m = λ / r
LAMBDA.deref().div(R.clone()));
pub const M_DASH: Lazy<BigUint> = Lazy::new(||
// m' = m/d
M.clone().div(D.clone()));

// e = 6X + 2;
pub static E: Lazy<BigUint> = Lazy::new(|| {
    let x6 = BigUint::from_i8(6).unwrap().mul(X.clone());
    // 6x + 2
    x6 + BigUint::from_i8(2).unwrap()
});

// optimal lambda in miller loop, lambda
pub const LAMBDA: Lazy<BigUint> = Lazy::new(|| {
    // lambdax = 6X + 2 + p - p^2 + p^3
    let p = MODULUS.clone();
    let p_pow2 = p.clone().pow(2_u32);
    let p_pow3 = p.clone().pow(3_u32);
    let x6 = BigUint::from_i8(6).unwrap().mul(X.deref());
    let two = BigUint::from_i8(2).unwrap();
    p_pow3.sub(p_pow2).add(p).add(two).add(x6)
});

pub static BETA: Lazy<Fq2> = Lazy::new(|| Fq2::new(Fq::from(9), Fq::ONE));

pub static BETA_PI_1: Lazy<Vec<Fq2>> = Lazy::new(|| {
    let mut res = vec![];
    for i in 1..6 {
        // exp = i * ((module-1)/6)
        let mut t = params::MODULUS.clone();
        t.sub_assign(BigUint::one());
        t.div_assign(BigUint::from_i32(6).unwrap());
        t.mul_assign(BigUint::from_i32(i).unwrap());
        let exp = t;
        let pi = BETA.pow(exp.to_u64_digits());
        res.push(pi);
    }
    res
});

pub static BETA_PI_2: Lazy<Vec<Fq2>> = Lazy::new(|| {
    let mut res = vec![];
    for i in 1..6 {
        // exp = i * ((module^2 -1)/6)
        let mut t = params::MODULUS.clone();
        t = t.pow(2_u32);
        t.sub_assign(BigUint::one());
        t.div_assign(BigUint::from_i32(6).unwrap());
        t.mul_assign(BigUint::from_i32(i).unwrap());
        let exp = t;
        let pi = BETA.pow(exp.to_u64_digits());
        res.push(pi);
    }
    res
});

pub static BETA_PI_3: Lazy<Vec<Fq2>> = Lazy::new(|| {
    let mut res = vec![];
    for i in 1..6 {
        // exp = i * ((module^3 -1)/6)
        let mut t = MODULUS.clone();
        t = t.pow(3_u32);
        t.sub_assign(BigUint::one());
        t.div_assign(BigUint::from_i32(6).unwrap());
        t.mul_assign(BigUint::from_i32(i).unwrap());
        let exp = t;
        let pi = BETA.pow(exp.to_u64_digits());
        res.push(pi);
    }
    res
});

pub fn tx(x: BigUint) -> BigUint {
    let p1 = BigUint::from_i8(6).unwrap();
    p1 * x.pow(2_u32) + BigUint::one()
}

pub fn hx(x: BigUint) -> BigUint {
    let p1 = BigUint::from_i8(18).unwrap();
    let numerator = p1 * x.pow(12_u32) - BigUint::one();
    let denominator = R.clone();
    numerator / denominator
}

pub fn mx(x: &BigUint) -> BigUint {
    let p1 = BigUint::from_i8(36).unwrap();
    let p2 = BigUint::from_i8(18).unwrap();
    p1.clone() * x.pow(4_u32) + p1 * x.pow(3_u32) + p2 * x.pow(2_u32) + BigUint::one()
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::utils::biguint_to_naf;

    use ark_ec::bn::BnConfig;
    use ark_ff::PrimeField;
    use std::str::FromStr;

    #[test]
    fn test_equivalently() {
        //  λ = 3rm′
        let actual = M_DASH
            .clone()
            .mul(R.clone())
            .mul(BigUint::from_i8(3).unwrap());
        assert_eq!(actual, LAMBDA.clone());

        // r_inv * r % h == 1
        let actual = R.clone().mul(R_INV.clone()) % &H.clone();
        assert_eq!(actual, BigUint::one());
        let expect = M.clone().mul(M_DASH.clone()) % &H.clone();
        assert_eq!(actual, BigUint::one());
    }

    #[test]
    fn test_lambda() {
        let actual = LAMBDA.clone();
        let expect = BigUint::from_str(
            "10486551571378427818905133077457505975146652579011797175399169355881771981095211883813744499745558409789005132135496770941292989421431235276221147148858384772096778432243207188878598198850276842458913349817007302752534892127325269"
        ).unwrap();

        assert_eq!(actual, expect);
    }

    #[test]
    fn test_modulus() {
        let actual = BigUint::from(Fq::MODULUS);
        let expect = BigUint::from_str_radix(MODULUS_STR, 16).unwrap();
        assert_eq!(actual, expect);

        let actual = BigUint::from_str(
            "21888242871839275222246405745257275088696311157297823662689037894645226208583",
        )
        .unwrap();
        assert_eq!(expect, actual);
    }

    #[test]
    fn test_e() {
        let actual = E.clone();
        let expect = BigUint::from_str("29793968203157093288").unwrap();
        assert_eq!(actual, expect);
    }

    #[test]
    fn test_ATE_LOOP_COUNT() {
        let digtals_naf = biguint_to_naf(E.clone());
        // NOTE:
        // Even if e_naf is different with ark_bn254::Config::ATE_LOOP_COUNT,
        // they play the same role in pairing_verifier.
        // assert_eq!(ark_bn254::Config::ATE_LOOP_COUNT, digtals_naf);
        let ATE_LOOP_COUNT_len = ark_bn254::Config::ATE_LOOP_COUNT.len();
        println!("ATE_LOOP_COUNT len: {:?}", ATE_LOOP_COUNT_len);
    }

    #[test]
    fn test_beta() {
        println!("Beta: {:?}", BETA.clone());

        println!("BETA_PI_1_2: {:?}", BETA_PI_1[1].to_string());
        println!("BETA_PI_1_3: {:?}", BETA_PI_1[2].to_string());
        println!("BETA_PI_2_2: {:?}", BETA_PI_2[1].to_string());
    }
}
