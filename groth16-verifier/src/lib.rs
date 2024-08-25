//! Normally Groth16 verifier did following:
//!  0. setup and generate Groth16's proof
//!  1. prepare the inputs: msm with vk and public_inputs
//!  2. check pairing with 4 pairings
//!      2.1 compute f with miller_loop
//!      2.2 compute final_f with `final_exponentiation`
//!
//! The pairing is costly. [`Prove on Pairing`] has figured out a solution: prove and verify for pairing.
//! By precompute miller lines and avoid the `final_exponentiation`, we can reduce the cost of paring.
//!
//! We can import the optimize of pairing to Groth16 verifier by recursive proof.
//! So the new Groth16 verifier does as following:
//!  0. setup and generate Groth16's proof
//!  1. prepare the inputs: msm with vk and public_inputs
//!  2. check pairing with 4 pairings
//!      2.1 setup for proving pairing: precompute lines and find_c.
//!      2.2 generate pairing proof
//!      2.2 verify the pairing proof.
//!
//!
//! Note: Only support Bn254 for now.
pub mod dummy_circuit;
#[cfg(test)]
mod test;

use ark_bn254::{Bn254, Fr, G1Affine, G1Projective};
use ark_ec::bn::G2Prepared;
use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, CurveGroup, VariableBaseMSM};
use ark_groth16::{PreparedVerifyingKey, Proof};
use ark_relations::r1cs::Result as R1CSResult;
use on_proving_pairings::prover::PairingProver;
use on_proving_pairings::setup::PairingPVKey;
use on_proving_pairings::verifier::PairingVerifier;
use std::ops::Neg;

pub struct Groth16Verifier;

impl Groth16Verifier {
    pub fn verify_proof_with_c_wi(
        pvk: &PreparedVerifyingKey<Bn254>,
        proof: &Proof<Bn254>,
        public_inputs: &[Fr],
    ) -> R1CSResult<bool> {
        let prepared_inputs = Self::prepare_inputs(pvk, public_inputs)?;
        Self::verify_proof_with_prepared_inputs(pvk, proof, &prepared_inputs)
    }

    // Porting from `ark_groth16::Groth16::prepare_inputs`
    pub fn prepare_inputs(
        pvk: &PreparedVerifyingKey<Bn254>,
        public_inputs: &[Fr],
    ) -> R1CSResult<G1Projective> {
        assert_eq!(public_inputs.len() + 1, pvk.vk.gamma_abc_g1.len());
        // g_ic = pvk.vk.gamma_abc_g1[0] + msm(pvk.vk.gamma_abc_g1[1..], public_inputs)
        let g_ic = pvk.vk.gamma_abc_g1[0].into_group();

        let g_ic = g_ic + G1Projective::msm(&pvk.vk.gamma_abc_g1[1..], &public_inputs).unwrap();

        Ok(g_ic)
    }

    // Verifier by applying with new paper:
    //
    // Porting from `ark_groth16::Groth16::verify_proof_with_prepared_inputs`
    pub fn verify_proof_with_prepared_inputs(
        pvk: &PreparedVerifyingKey<Bn254>,
        proof: &Proof<Bn254>,
        prepared_inputs: &G1Projective,
    ) -> R1CSResult<bool> {
        let beta_prepared: G2Prepared<ark_bn254::Config> = (pvk.vk.beta_g2.clone().neg()).into();
        let gamma_g2_neg_pc = pvk.gamma_g2_neg_pc.clone();
        let delta_g2_neg_pc = pvk.delta_g2_neg_pc.clone();
        let sum_ai_abc_gamma = prepared_inputs.into_affine();

        // Pi
        let a = vec![
            sum_ai_abc_gamma.into(),
            proof.c.into(),
            pvk.vk.alpha_g1.into(),
            <G1Affine as Into<<Bn254 as Pairing>::G1Prepared>>::into(proof.a),
        ];
        // Qi - proof.a(Q4), is non-fixed.
        let b = vec![
            pvk.gamma_g2_neg_pc.clone(),
            pvk.delta_g2_neg_pc.clone(),
            beta_prepared.clone(),
            proof.b.into(),
        ];

        // Prove and verify pairing.
        // finding_c
        let pairing_pvk = PairingPVKey::setup(a.clone(), b.clone());

        // eval_points: [P1,P2,P3]
        let eval_points = vec![sum_ai_abc_gamma, proof.c, pvk.vk.alpha_g1];
        // precompute lines: [Q1,Q2,Q3]
        let q_prepared_lines = b[0..3].to_vec();

        let final_f = PairingProver::prove_quad_pairing(
            eval_points,
            &q_prepared_lines,
            proof.a,
            proof.b,
            &pairing_pvk,
        );

        // verify
        PairingVerifier::verify(&pairing_pvk, final_f);
        Ok(true)
    }
}
