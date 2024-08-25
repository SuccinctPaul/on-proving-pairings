use num_bigint::BigUint;
use num_traits::ToPrimitive;

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

#[cfg(test)]
mod test {
    use super::*;
    use crate::params;

    use ark_bn254::G2Affine;
    use ark_ec::bn::g2::mul_by_char;
    use ark_std::UniformRand;
    use std::ops::Mul;
    use std::ops::{Deref, Neg};

    #[test]
    fn test_biguint_naf() {
        let mut expect = to_naf(29793968203157093288);
        expect.reverse();
        expect.remove(0);
        println!("res: {:?}", expect);

        println!("E :{:?}", params::E.deref());
        let mut actual = crate::utils::biguint_to_naf(params::E.clone());
        actual.reverse();
        actual.remove(0);
        println!("res: {:?}", actual);

        assert_eq!(expect, actual);
    }

    // #[test]
    // fn test_home_projective_and_projective(){
    //     // TODO Do more test for this.
    //     let mut rng = test_rng();
    //
    //     let p = G1Affine::rand(&mut rng);
    //     let t = G2Affine::rand(&mut rng).into_group();
    //     let q = G2Affine::rand(&mut rng);
    //
    //     let mut expect = G2HomProjective {
    //         x: t.x,
    //         y: t.y,
    //         z: t.z,
    //     };
    //     expect.add_in_place(&q);
    // }

    #[test]
    fn test_mul_by_char() {
        let rng = &mut ark_std::test_rng();
        let Q4 = G2Affine::rand(rng);

        // ==== a. Compute phi_Q
        let actual = mul_by_char::<ark_bn254::Config>(Q4.clone());

        // 1. one-time frobenius map to compute phi_Q
        // 2.1 Qx.conjugate * beta^{2 * (p - 1) / 6}
        let mut Q4_x = Q4.x.clone();
        Q4_x.conjugate_in_place();
        Q4_x = Q4_x.mul(&params::BETA_PI_1[1]);

        // 2.2 Qy.conjugate * beta^{3 * (p - 1) / 6}
        let mut Q4_y = Q4.y.clone();
        Q4_y.conjugate_in_place();
        Q4_y = Q4_y.mul(&params::BETA_PI_1[2]);
        let phi_Q = G2Affine::new(Q4_x, Q4_y);

        assert_eq!(actual, phi_Q);

        // ==== B. Compute phi_Q_2
        let mut actual_2 = mul_by_char::<ark_bn254::Config>(phi_Q.clone());
        actual_2.y = actual_2.y.neg();

        // 2. one-time frobenius map to compute phi_Q
        // 2.1  Qx * beta^{2 * (p^2 - 1) / 6}
        let mut Q4_x = Q4.x.clone();
        Q4_x = Q4_x.mul(&params::BETA_PI_2[1]);

        // 2.2 Qy
        // let Q4_y = Q4.y.clone().neg();
        let Q4_y = Q4.y.clone();

        // 2.3(non-fixed) add line with T4 and phi(Q)
        let phi_Q_2 = G2Affine::new(Q4_x, Q4_y);

        assert_eq!(phi_Q_2.y, Q4.y);
        assert_eq!(actual_2.y, Q4.y);
        assert_eq!(actual_2, phi_Q_2);
    }
}
