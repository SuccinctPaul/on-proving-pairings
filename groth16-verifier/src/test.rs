use crate::Groth16Verifier;

use crate::dummy_circuit::gen_groth16_dummy_circuit_proof;
use ark_bn254::Bn254;
use ark_groth16::Groth16;

#[test]
fn test_groth16_verifier() {
    type E = Bn254;

    let k = 6;

    // 1. gen proof
    let (proof, pvk, pi) = gen_groth16_dummy_circuit_proof::<E>(k);

    // 2. verify with native verifier
    assert!(
        Groth16::<E>::verify_proof(&pvk, &proof, &pi).unwrap(),
        "native verifier can't pass"
    );

    // 3. verifier with new one
    assert!(Groth16Verifier::verify_proof_with_c_wi(&pvk, &proof, &pi).unwrap());
}
