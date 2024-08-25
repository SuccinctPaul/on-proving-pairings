use crate::lambda_residues::LambdaResidues;
use ark_bn254::{Bn254, Fq12};
use ark_ec::bn::{G1Prepared, G2Prepared};
use ark_ec::pairing::Pairing;
use ark_ff::Field;

pub struct PairingPVKey {
    // f, then base line
    pub f: Fq12,
    // c
    pub c: Fq12,
    pub wi: Fq12,
    pub c_inv: Fq12,
}

impl PairingPVKey {
    pub fn setup(
        a: impl IntoIterator<Item = impl Into<G1Prepared<ark_bn254::Config>>>,
        b: impl IntoIterator<Item = impl Into<G2Prepared<ark_bn254::Config>>>,
    ) -> Self {
        // compute f, then base line
        let qap = Bn254::multi_miller_loop(a, b);
        let f = qap.0;

        // finding_c
        let witness = LambdaResidues::finding_c(f);
        let c_inv = witness.c.inverse().unwrap();

        Self {
            f,
            c: witness.c,
            wi: witness.wi,
            c_inv,
        }
    }
}
