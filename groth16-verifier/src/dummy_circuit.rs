use ark_crypto_primitives::snark::{CircuitSpecificSetupSNARK, SNARK};
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use ark_groth16::{prepare_verifying_key, Groth16};
use ark_relations::lc;
use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError};
use ark_std::{test_rng, UniformRand};
use rand::{RngCore, SeedableRng};

#[derive(Copy)]
struct DummyCircuit<F: PrimeField> {
    pub a: Option<F>,
    pub b: Option<F>,
    pub num_variables: usize,
    pub num_constraints: usize,
}

impl<F: PrimeField> Clone for DummyCircuit<F> {
    fn clone(&self) -> Self {
        DummyCircuit {
            a: self.a.clone(),
            b: self.b.clone(),
            num_variables: self.num_variables.clone(),
            num_constraints: self.num_constraints.clone(),
        }
    }
}

impl<F: PrimeField> ConstraintSynthesizer<F> for DummyCircuit<F> {
    fn generate_constraints(self, cs: ConstraintSystemRef<F>) -> Result<(), SynthesisError> {
        let a = cs.new_witness_variable(|| self.a.ok_or(SynthesisError::AssignmentMissing))?;
        let b = cs.new_witness_variable(|| self.b.ok_or(SynthesisError::AssignmentMissing))?;
        let c = cs.new_input_variable(|| {
            let a = self.a.ok_or(SynthesisError::AssignmentMissing)?;
            let b = self.b.ok_or(SynthesisError::AssignmentMissing)?;

            Ok(a * b)
        })?;

        for _ in 0..(self.num_variables - 3) {
            let _ = cs.new_witness_variable(|| self.a.ok_or(SynthesisError::AssignmentMissing))?;
        }

        for _ in 0..self.num_constraints - 1 {
            cs.enforce_constraint(lc!() + a, lc!() + b, lc!() + c)?;
        }

        cs.enforce_constraint(lc!(), lc!(), lc!())?;

        Ok(())
    }
}

// Gen a groth16-verifier proof to verifier
//
// return: (proof, pvk, public_inputs)
pub fn gen_groth16_dummy_circuit_proof<E: Pairing>(
    k: usize,
) -> (
    ark_groth16::Proof<E>,
    ark_groth16::PreparedVerifyingKey<E>,
    Vec<E::ScalarField>,
) {
    let mut rng = ark_std::rand::rngs::StdRng::seed_from_u64(test_rng().next_u64());

    let circuit = DummyCircuit::<E::ScalarField> {
        a: Some(E::ScalarField::rand(&mut rng)),
        b: Some(E::ScalarField::rand(&mut rng)),
        num_variables: 10,
        num_constraints: 1 << k,
    };

    let (pk, vk) = Groth16::<E>::setup(circuit, &mut rng).unwrap();
    let pvk = prepare_verifying_key::<E>(&vk);

    let c = circuit.a.unwrap() * circuit.b.unwrap();

    let proof = Groth16::<E>::prove(&pk, circuit, &mut rng).unwrap();

    // public inputs
    let pi = vec![c];

    (proof, pvk, pi)
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bn254::Bn254;

    use ark_groth16::Groth16;

    #[test]
    fn test_gen_groth16_proof() {
        type E = Bn254;

        let k = 6;
        let (proof, pvk, pi) = gen_groth16_dummy_circuit_proof::<E>(k);
        assert!(Groth16::<E>::verify_proof(&pvk, &proof, &pi).unwrap());
    }
}
