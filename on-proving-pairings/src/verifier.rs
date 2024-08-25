use crate::params;
use crate::setup::PairingPVKey;
use ark_bn254::Fq12;
use ark_ff::Field;

pub struct PairingVerifier;

impl PairingVerifier {
    pub fn verify(pvk: &PairingPVKey, final_f: Fq12) {
        let p_pow3 = params::MODULUS.pow(3_u32);
        let lambda = params::LAMBDA.clone();
        let (exp, sign) = if lambda > p_pow3 {
            (lambda - p_pow3, true)
        } else {
            (p_pow3 - lambda, false)
        };

        let hint = if sign {
            pvk.f * pvk.wi * (pvk.c_inv.pow(exp.to_u64_digits()))
        } else {
            pvk.f * pvk.wi * (pvk.c_inv.pow(exp.to_u64_digits()).inverse().unwrap())
        };

        let p_pow3 = params::MODULUS.pow(3_u32);
        assert_eq!(hint, pvk.c.pow(p_pow3.to_u64_digits()), "hint is wrong");
        assert_eq!(final_f, hint, "final_f not equal hint");
    }
}
