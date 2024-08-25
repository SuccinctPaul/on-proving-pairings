use crate::params;
use crate::params::{E, LAMBDA, MODULUS};
use crate::utils::biguint_to_naf;
use ark_bn254::{Fq, Fq12, Fq2, Fq6, G1Affine, G2Affine, G2Projective};
use ark_ec::{AffineRepr, CurveGroup, Group};
use ark_ff::Field;
use num_bigint::BigUint;
use num_traits::Zero;
use std::ops::{Add, Div, Mul, Neg};

// stands for (alpha, bias)
pub type MillerLinesRes = (Fq2, Fq2);

// A significant part of Miller function evaluation consists of line constructions and evaluations
pub struct MillerLines;
// TODO Replace it later.
// {
//     alpha: Fq2,
//     bias: Fq2,
// }

impl MillerLines {
    pub fn mul_line_base(r: Fq12, a: Fq2, b: Fq2, c: Fq2) -> Fq12 {
        let fq6_y = Fq6::new(b, a, Fq2::ZERO);
        let fq6_x = Fq6::new(c, Fq2::ZERO, Fq2::ZERO);
        let fl = Fq12::new(fq6_x, fq6_y);
        r.mul(fl)
    }

    // we use affine coordinate to verify the line evaluation
    // (-b) + y_P * w^3 + (-alpha * x_P) * w^2 where w \in Fp12
    pub fn line_evaluation(alpha: Fq2, bias: Fq2, point: G1Affine) -> (Fq2, Fq2, Fq2) {
        let mut neg_alpha = alpha.neg();
        neg_alpha.mul_assign_by_basefield(&point.x);
        (bias.neg(), neg_alpha, Fq2::new(point.y, Fq::ZERO))
    }

    fn line_double(point: &G2Affine) -> MillerLinesRes {
        let (x, y) = (point.x, point.y);

        // slope: alpha = 3 * x ^ 2 / (2 * y)
        let alpha = x.square().mul(Fq2::from(3)).div(y.mul(Fq2::from(2)));
        // bias = y - alpha * x
        let bias = y - alpha * x;

        (alpha, bias)
    }

    // NOTE: point can't equal with other.
    fn line_add(point: &G2Affine, other: &G2Affine) -> MillerLinesRes {
        let (x1, y1) = (point.x, point.y);
        let (x2, y2) = (other.x, other.y);

        // slope: alpha = (y2 - y1) / (x2 - x1)
        let alpha = (y2 - y1) / (x2 - x1);
        // bias = y1 - alpha * x1
        let bias = y1 - alpha * x1;
        (alpha, bias)
    }

    // Indexer for fixed point Q
    // cache line line_function for [6x + 2 + p - p^2 + p^3]Q
    // Input:
    //      Q: The Fixed Point
    // Ref: Algorithm 7 of [On Proving Pairings](https://eprint.iacr.org/2024/640.pdf)
    pub fn precompute_lines(Q: G2Projective, e: BigUint, lamb: BigUint) -> Vec<MillerLinesRes> {
        let mut naf_digits = biguint_to_naf(e.clone());
        naf_digits.reverse();
        naf_digits.remove(0);

        let mut line_vec = vec![];

        let mut T = Q.clone();

        // 1. double-add part, 6x + 2
        naf_digits.iter().enumerate().for_each(|(i, digit)| {
            let double_res = MillerLines::line_double(&T.into_affine());
            line_vec.push(double_res);

            T = T.double();
            let digtil_pow2 = digit * digit;
            if digtil_pow2 == 1 {
                let qt = if 1 == *digit {
                    Q.clone()
                } else {
                    Q.clone().neg()
                };

                let qt_double_res = Self::line_add(&T.into_affine(), &qt.into_affine());
                T = T.add(qt);
                line_vec.push(qt_double_res);
            }
        });
        assert_eq!(T, Q.into_affine().mul_bigint(e.to_u64_digits()));

        // 2. frobenius map part, p - p^2 + p^3
        // 2.1 Q1 = pi(Q)
        // x = x' * beta^(2 * (p - 1) / 6)
        // y = y' * beta^(3 * (p - 1) / 6))
        let (mut x, mut y) = (Q.x.clone(), Q.y.clone());

        let pi_1_Q = G2Projective::new(
            x.conjugate_in_place().mul(params::BETA_PI_1[1]),
            y.conjugate_in_place().mul(params::BETA_PI_1[2]),
            Fq2::ONE,
        );
        assert_eq!(pi_1_Q, Q.into_affine().mul_bigint(MODULUS.to_u64_digits()));

        // 2.2. Q2 = pi2(Q)
        // x = x * beta * (2 * (p^2 - 1) / 6)
        // y = y * beta * (3 * (p^2 - 1) / 6) = -y
        let (x, y) = (Q.x, Q.y);
        let pi_2_Q = G2Projective::new(
            x.mul(params::BETA_PI_2[1]),
            y.mul(params::BETA_PI_2[2]),
            Fq2::ONE,
        );
        assert_eq!(
            pi_2_Q,
            Q.into_affine().mul_bigint(MODULUS.pow(2).to_u64_digits())
        );

        // 2.3. Q3 = pi3(Q)
        // x = x' * beta * (2 * (p^3 - 1) / 6)
        // y = y' * beta * (3 * (p^3 - 1) / 6)
        let (mut x, mut y) = (Q.x.clone(), Q.y.clone());

        let pi_3_Q = G2Projective::new(
            x.conjugate_in_place().mul(params::BETA_PI_3[1]),
            y.conjugate_in_place().mul(params::BETA_PI_3[2]),
            Fq2::ONE,
        );
        assert!(pi_3_Q.into_affine().is_on_curve());
        assert_eq!(
            pi_3_Q,
            Q.into_affine().mul_bigint(MODULUS.pow(3).to_u64_digits())
        );

        let line_pi_1 = Self::line_add(&T.into_affine(), &pi_1_Q.into_affine());
        T = T.add(pi_1_Q);
        line_vec.push(line_pi_1);
        assert_eq!(
            T,
            Q.into_affine()
                .mul_bigint(E.clone().add(MODULUS.clone()).to_u64_digits())
        );

        let neg_pi_2_Q = pi_2_Q.neg();

        let line_pi_2 = Self::line_add(&T.into_affine(), &neg_pi_2_Q.into_affine());
        T = T.add(neg_pi_2_Q);

        line_vec.push(line_pi_2);

        // TODO: Move this outersides
        //    k = e + px(x) - px(x) ** 2
        //     assert(T == Q.scalar_mul(k if k > 0 else rx(x) - (-k % rx(x))))
        // k = p - p^2 + e
        // let k = MODULUS.clone().add(E.clone()).sub(MODULUS.clone().pow(2)) ;
        // let k = k + e;
        // assert_eq!(T, Q.mul_bigint(if k.gt(BigUint::ZERO){
        //     k
        // }else {
        // TODO
        //     // rx(x) - (-k % rx(x))
        // } ));

        let line_pi_3 = (Fq2::ZERO, T.x.mul(T.z.inverse().unwrap().square()));
        T = T.add(pi_3_Q);
        line_vec.push(line_pi_3);

        assert!(line_pi_3.0.is_zero());
        assert!(T.into_affine().is_on_curve());
        assert_eq!(T, Q.mul_bigint(LAMBDA.to_u64_digits()));

        line_vec
    }
}

#[cfg(test)]
mod test {
    use ark_bn254::{Fq, Fq2};
    use ark_ec::{AffineRepr, CurveGroup};
    use ark_std::{end_timer, start_timer};
    use num_bigint::BigUint;
    use num_traits::{FromPrimitive, One};
    use std::str::FromStr;

    use super::*;
    use crate::{
        dev::{g1, g2},
        params::E,
    };

    #[test]
    fn test_line_double() {
        println!("Start");
        let Q1 = G2Projective::from(g2);

        println!("Q1: {:?}", Q1);
        let line = MillerLines::line_double(&Q1.into_affine());

        println!("line: {:?}", line.0.to_string());
        println!("line: {:?}", line.1.to_string());

        assert_eq!(
            line.0,
            Fq2::new(
                Fq::from_str(
                    "3011624467519903477642955996195800747228099871883200115829824046966119461903"
                )
                .unwrap(),
                Fq::from_str(
                    "8108057915722081371595498492026831474963176496200081396505507133092048631529"
                )
                .unwrap(),
            )
        );

        assert_eq!(
            line.1,
            Fq2::new(
                Fq::from_str(
                    "5093819629657497304798450698428222859493999469989845044032174107566603880340"
                )
                .unwrap(),
                Fq::from_str(
                    "1206257300307291756486735540503199561073344972686232605282019777893274262471"
                )
                .unwrap(),
            )
        );
    }

    #[test]
    fn test_line_add() {
        println!("Start");
        let Q1 = G2Projective::from(g2);
        let Q2 = G2Projective::from(g2).double();

        println!("Q1: {:?}", Q1.into_affine());
        println!("Q2: {:?}", Q2.into_affine());
        let line = MillerLines::line_add(&Q1.into_affine(), &Q2.into_affine());
        // println!("alpha: {:?}", line.0.to_string());
        println!("bias: {:?}", line.1.to_string());
        assert_eq!(
            line.0,
            Fq2::new(
                Fq::from_str(
                    "17322200159135559786593682602367196430903128621957194838174337471042229537489"
                )
                .unwrap(),
                Fq::from_str(
                    "13357879408196058558132338895939622053982269645137478732965497231920858080670"
                )
                .unwrap(),
            )
        );

        assert_eq!(
            line.1,
            Fq2::new(
                Fq::from_str(
                    "9837371316849015087103862504041462238679211850457142473950276246745360727711"
                )
                .unwrap(),
                Fq::from_str(
                    "20187872875180194894625542094089895199219468421059757009464967542497725908996"
                )
                .unwrap(),
            )
        );
    }

    #[test]
    fn test_line_precomputation() {
        // assume we want to prove e(P1, Q1) = e(P2, Q2), namely e(P1, Q1) * e(P2, -Q2) = 1
        // fixed point Q in G2, public known to verifier
        let P1 = g1.mul_bigint(BigUint::from_i8(3).unwrap().to_u64_digits());
        let P2 = g1.mul_bigint(BigUint::one().to_u64_digits());

        let Q1 = g2.mul_bigint(BigUint::one().to_u64_digits());
        let Q2 = g2.mul_bigint(BigUint::from_i8(3).unwrap().to_u64_digits());

        // indexer (oracle) for fixed point Q
        println!("Q1: {:?}", Q1);
        println!("e: {:?}", E.clone());
        println!("LAMBDA: {:?}", LAMBDA.clone());
        let start = start_timer!(|| "start compute line precomputation");
        let L1 = MillerLines::precompute_lines(Q1, E.clone(), LAMBDA.clone());
        end_timer!(start);
    }
}
