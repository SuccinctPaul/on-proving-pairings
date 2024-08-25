/// Generate proof of pairings.
use crate::params;
use crate::setup::PairingPVKey;
use ark_bn254::{Bn254, Fq12, G1Affine, G2Affine};
use ark_ec::bn::g2::{mul_by_char, G2HomProjective};
use ark_ec::bn::{BnConfig, G2Prepared};
use ark_ec::AffineRepr;
use ark_ff::Field;
use num_traits::One;
use std::ops::Neg;

// Prove pairing and generate proof: final_f.
//
// And can verify the proof by check:  c^lambda = final_f * wi, namely c_inv^lambda * final_f * wi = 1
pub struct PairingProver;

impl PairingProver {
    // Prove pairing with duality pairing.
    // todo

    // Prove quadratic pairing.
    // Will generate proof: final_f to prove `e(P1,Q1)*e(P2,Q2)*e(P3,Q3)*e(P4,Q4)=1`
    //
    // Note: here only supports to prove on `groth16-verifier's pairing`,
    //          where (Q1,Q2,Q3) are fixed, and Q4 is non-fixed.
    //
    // params:
    //  @eval_points: [P1,P2,P3]. which maps with fixed [Q1,Q2,Q3]
    //  @lines: precompute miller lines for Qi. Only support fixed Qi.
    //          Base line: [Q1,Q2,Q3]
    //  @P4: P4
    //  @Q4: non fixed base line.
    //  @c: c^lambda = f*w^i
    //  @c_inv: inverse of c
    //
    // verify c^lambda = f * wi, namely c_inv^lambda * f * wi = 1
    pub fn prove_quad_pairing(
        eval_points: Vec<G1Affine>,
        lines: &Vec<G2Prepared<ark_bn254::Config>>,
        P4: G1Affine,
        Q4: G2Affine,
        pvk: &PairingPVKey,
    ) -> Fq12 {
        assert_eq!(eval_points.len(), 3, "Should contains 3 G1Affine: P1,P2,P3");
        assert_eq!(lines.len(), 3, "Only 3 precompute lines for Q1,Q2,Q3");
        assert_eq!(pvk.c * pvk.c_inv, Fq12::ONE, "Check if c·c^−1 = 1");

        let mut T4 = G2HomProjective::<ark_bn254::Config> {
            x: Q4.x,
            y: Q4.y,
            z: ark_bn254::Fq2::one(),
        };

        // constants: 1/2
        let two_inv = ark_bn254::Fq::one().double().inverse().unwrap();

        // 1. f = c_inv
        let mut f = pvk.c_inv;

        let mut lines_iters = lines
            .iter()
            .map(|item| item.ell_coeffs.iter())
            .collect::<Vec<_>>();

        // 2. miller loop part, 6x + 2
        for i in (1..ark_bn254::Config::ATE_LOOP_COUNT.len()).rev() {
            let bit = ark_bn254::Config::ATE_LOOP_COUNT[i - 1];

            // 2.1 double: f = f * f
            f = f.square();

            // 2.2 mul c
            //  f = f * c_inv, if digit == 1
            //  f = f * c, if digit == -1
            f = if 1 == bit {
                f * pvk.c_inv
            } else if bit == -1 {
                f * pvk.c
            } else if bit == 0 {
                f
            } else {
                panic!("bit is not in (-1,1), bit={bit}");
            };

            // 2.3 accumulate double lines (fixed and non-fixed)
            // 2.3.1(fixed) f = f * double_line_Q(P). fixed points: P1, P2, P3
            for (line_i, pi) in lines_iters.iter_mut().zip(eval_points.iter()) {
                let line_i_0 = line_i.next().unwrap();
                Bn254::ell(&mut f, line_i_0, pi);
            }

            // 2.3.2(non-fixed) double line with T4 (projective coordinates)
            let double_line = T4.double_in_place(&two_inv); // TODO: check if the param is 1/2

            // 2.3.3(non-fixed) evaluation double_line. non-fixed points: P4
            Bn254::ell(&mut f, &double_line, &P4);

            if bit == 1 || bit == -1 {
                // 2.4 accumulate add lines (fixed and non-fixed)
                // 2.4.1(fixed) f = f * add_line_eval. fixed points: P1, P2, P3
                for (line_i, pi) in lines_iters.iter_mut().zip(eval_points.iter()) {
                    let line_i_1 = line_i.next().unwrap();
                    Bn254::ell(&mut f, line_i_1, pi);
                }
                // 2.4.2(non-fixed) double line with T4 (projective coordinates)
                let add_line = if bit == 1 {
                    T4.add_in_place(&Q4)
                } else {
                    // }else if bit == -1 {
                    T4.add_in_place(&Q4.neg())
                };

                // 2.4.3(non-fixed) evaluation double_line. non-fixed points: P4
                Bn254::ell(&mut f, &add_line, &P4);
            }
        }

        // 3. f = f * c_inv^p * c^{p^2}
        f = f
            * pvk.c_inv.pow(params::MODULUS.to_u64_digits())
            * pvk.c.pow(params::MODULUS.pow(2).to_u64_digits());

        // 4. scale f: f = f * wi .
        f = f * pvk.wi;

        // 5 add lines (fixed and non-fixed)
        // 5.1(fixed) f = f * add_line_eval. fixed points: P1, P2, P3
        for (line_i, pi) in lines_iters.iter_mut().zip(eval_points.iter()) {
            let line_i_1 = line_i.next().unwrap();
            Bn254::ell(&mut f, line_i_1, pi);
        }
        // 5.2(non-fixed) one-time frobenius map to compute phi_Q
        //     compute phi(Q) with Q4
        let phi_Q = mul_by_char::<ark_bn254::Config>(Q4.clone());

        // 5.3(non-fixed) add line with phi_Q
        let add_line = T4.add_in_place(&phi_Q);

        // 5.4(non-fixed) evaluation add_lin. non-fixed points: P4
        Bn254::ell(&mut f, &add_line, &P4);

        // 6. add lines (fixed and non-fixed)
        // 6.1(fixed) f = f * add_line_eval. fixed points: P1, P2, P3
        for (line_i, pi) in lines_iters.iter_mut().zip(eval_points.iter()) {
            let line_i_1 = line_i.next().unwrap();
            Bn254::ell(&mut f, line_i_1, pi);
        }
        // 6.2 two-time frobenius map to compute phi_Q
        //     compute phi_Q_2 with phi_Q
        // mul_by_char: used to q's frob...map.
        let mut phi_Q_2 = mul_by_char::<ark_bn254::Config>(phi_Q.clone());
        phi_Q_2.y.neg_in_place();

        // 6.3 add line with phi_Q_2
        let add_line = T4.add_in_place(&phi_Q_2);

        // 6.4 evaluation add_lin. non-fixed points: P4
        Bn254::ell(&mut f, &add_line, &P4);

        // return final_f
        f
    }
}
