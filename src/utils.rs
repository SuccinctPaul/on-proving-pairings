use ark_ff::One;
use num_bigint::BigUint;
use num_traits::{FromPrimitive, ToPrimitive};

pub fn biguint_to_naf(num: BigUint) -> Vec<i8> {
    to_naf(num.to_i128().unwrap())
}

fn to_naf(mut x: i128) -> Vec<i8> {
    let mut z = vec![];
    while x > 0 {
        if x % 2 == 0 {
            z.push(0);
        } else {
            let zi: i8 = 2 - (x % 4) as i8;
            x -= zi as i128;
            z.push(zi)
        }

        x = x / 2
    }
    return z;
}

use crate::params;
use ark_bn254::{Fq12, Fq6};
use std::ops::Mul;

pub fn fq12_to_frobenius(mut q12: Fq12) -> Fq12 {
    let e2_z = q12.c1.c2.conjugate_in_place().mul(params::BETA_PI_1[4]);
    let e2_y = q12.c1.c1.conjugate_in_place().mul(params::BETA_PI_1[2]);
    let e2_x = q12.c1.c0.conjugate_in_place().mul(params::BETA_PI_1[0]);

    let e1_z = q12.c0.c2.conjugate_in_place().mul(params::BETA_PI_1[3]);
    let e1_y = q12.c0.c1.conjugate_in_place().mul(params::BETA_PI_1[1]);
    let e1_x = q12.c0.c0.conjugate_in_place().to_owned();

    Fq12::new(Fq6::new(e1_x, e1_y, e1_z), Fq6::new(e2_x, e2_y, e2_z))
}

pub fn fq12_to_frobenius_p2(q12: Fq12) -> Fq12 {
    let e2_z = q12.c1.c2 * params::BETA_PI_2[4];
    let e2_y = q12.c1.c1 * params::BETA_PI_2[2];
    let e2_x = q12.c1.c0 * params::BETA_PI_2[0];

    let e1_z = q12.c0.c2 * params::BETA_PI_2[3];
    let e1_y = q12.c0.c1 * params::BETA_PI_2[1];
    let e1_x = q12.c0.c0;

    Fq12::new(Fq6::new(e1_x, e1_y, e1_z), Fq6::new(e2_x, e2_y, e2_z))
}

pub fn fq12_to_frobenius_p3(mut q12: Fq12) -> Fq12 {
    let e2_z = q12.c1.c2.conjugate_in_place().mul(params::BETA_PI_3[4]);
    let e2_y = q12.c1.c1.conjugate_in_place().mul(params::BETA_PI_3[2]);
    let e2_x = q12.c1.c0.conjugate_in_place().mul(params::BETA_PI_3[0]);
    let e1_z = q12.c0.c2.conjugate_in_place().mul(params::BETA_PI_3[3]);
    let e1_y = q12.c0.c1.conjugate_in_place().mul(params::BETA_PI_3[1]);
    let e1_x = q12.c0.c0.conjugate_in_place().to_owned();

    Fq12::new(Fq6::new(e1_x, e1_y, e1_z), Fq6::new(e2_x, e2_y, e2_z))
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::params;
    use std::ops::Deref;

    #[test]
    fn test_biguint_naf() {
        let mut expect = to_naf(29793968203157093288);
        expect.reverse();
        expect.remove(0);
        println!("res: {:?}", expect);

        println!("E :{:?}", params::E.deref());
        let mut actual = biguint_to_naf(params::E.clone());
        actual.reverse();
        actual.remove(0);
        println!("res: {:?}", actual);

        assert_eq!(expect, actual);
    }

    #[test]
    fn test_beta_pi() {
        for x in params::BETA_PI_1.deref() {
            println!("beta_pi_1: {:?}", x.to_string());
        }
        println!("");
        for x in params::BETA_PI_1.deref() {
            println!("beta_pi_2: {:?}", x.to_string());
        }
        println!("");
        for x in params::BETA_PI_3.deref() {
            println!("beta_pi_3: {:?}", x.to_string());
        }
    }
}
