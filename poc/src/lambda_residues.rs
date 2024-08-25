use crate::params;
use crate::params::MODULUS;
use ark_bn254::Fq12;
use ark_ff::Field;
use ark_std::UniformRand;
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive};
use rand_chacha::rand_core::SeedableRng;
use rand_chacha::ChaCha20Rng;
use std::ops::Deref;

// c and wi,
// satisfying c^lambda = f * wi
pub struct LambdaResidues {
    pub c: Fq12,
    pub wi: Fq12,
}

impl LambdaResidues {
    // Computing λ residues over BN curve
    // Input:
    //      f: output of a Miller loop.
    //          It's always be r-th and m′-th residue, but it might not be a cubic residue.
    // Output:
    //      c and wi,
    //      satisfying c^lambda = f * wi
    //
    // Ref: Algorithm 5 of [On Proving Pairings](https://eprint.iacr.org/2024/640.pdf)
    pub fn finding_c(f: Fq12) -> Self {
        let s = 3_u32;
        let exp = MODULUS.pow(12_u32) - 1_u32;
        let h = &exp / params::R.deref();
        let t = &exp / 3_u32.pow(s);
        let k = (&t + 1_u32) / 3_u32;
        let m = params::LAMBDA.deref() / params::R.deref();
        let d = 3_u32;
        let mm = &m / d;

        let mut prng = ChaCha20Rng::seed_from_u64(0);
        let cofactor_cubic = 3_u32.pow(s - 1) * &t;

        // Find C. See more: 4.3.2 Finding c
        // make f is r-th residue, but it's not cubic residue
        assert_eq!(f.pow(h.to_u64_digits()), ark_bn254::Fq12::ONE);
        assert_ne!(f.pow(cofactor_cubic.to_u64_digits()), ark_bn254::Fq12::ONE);

        // sample a proper scalar w which is cubic non-residue
        let w = {
            let (mut w, mut z) = (ark_bn254::Fq12::ONE, ark_bn254::Fq12::ONE);
            while w == ark_bn254::Fq12::ONE {
                // choose z which is 3-th non-residue
                let mut legendre = ark_bn254::Fq12::ONE;
                while legendre == ark_bn254::Fq12::ONE {
                    z = ark_bn254::Fq12::rand(&mut prng);
                    legendre = z.pow(cofactor_cubic.to_u64_digits());
                }
                // obtain w which is t-th power of z
                w = z.pow(t.to_u64_digits());
            }
            w
        };
        // make sure 27-th root w, is 3-th non-residue and r-th residue
        assert_ne!(w.pow(cofactor_cubic.to_u64_digits()), ark_bn254::Fq12::ONE);
        assert_eq!(w.pow(h.to_u64_digits()), ark_bn254::Fq12::ONE);

        // just two option, w and w^2, since w^3 must be cubic residue, leading f*w^3 must not be cubic residue
        let mut wi = w;
        if (f * wi).pow(cofactor_cubic.to_u64_digits()) != ark_bn254::Fq12::ONE {
            assert_eq!(
                (f * w * w).pow(cofactor_cubic.to_u64_digits()),
                ark_bn254::Fq12::ONE
            );
            wi = w * w;
        }
        assert_eq!(wi.pow(h.to_u64_digits()), ark_bn254::Fq12::ONE);

        assert_eq!(params::LAMBDA.clone(), &d * &mm * params::R.deref());
        // f1 is scaled f
        let f1 = f * wi;

        // r-th root of f1, say f2
        let r_inv = params::R.deref().modinv(&h).unwrap();
        assert_ne!(r_inv, BigUint::one());
        let f2 = f1.pow(r_inv.to_u64_digits());
        assert_ne!(f2, ark_bn254::Fq12::ONE);

        // m'-th root of f, say f3
        let mm_inv = mm.modinv(&(params::R.deref() * h)).unwrap();
        assert_ne!(mm_inv, BigUint::one());
        let f3 = f2.pow(mm_inv.to_u64_digits());
        assert_eq!(f3.pow(cofactor_cubic.to_u64_digits()), ark_bn254::Fq12::ONE);
        assert_ne!(f3, ark_bn254::Fq12::ONE);

        // d-th (cubic) root, say c
        let c = Self::tonelli_shanks_cubic(f3, w, s, t, k);
        assert_ne!(c, ark_bn254::Fq12::ONE);
        assert_eq!(c.pow(params::LAMBDA.clone().to_u64_digits()), f * wi);

        Self { c, wi }
    }

    // Tonelli-Shanks for cube roots
    // Inputs:
    //         a: Fp12 which is cubic residue
    //         c: random Fp12 which is cubic non-residue
    //         s: satisfying p^12 - 1 = 3^s * t
    //         t: satisfying p^12 - 1 = 3^s * t
    //         k: k = (t + 1) // 3
    //
    // Ref: Table 3 from https://eprint.iacr.org/2009/457.pdf
    // Ref: Algorithm 4 of [On Proving Pairings](https://eprint.iacr.org/2024/640.pdf)
    pub fn tonelli_shanks_cubic(a: Fq12, c: Fq12, s: u32, t: BigUint, k: BigUint) -> Fq12 {
        // r = a^t
        let mut r = a.pow(t.to_u64_digits());
        let e = 3_u32.pow(s - 1);
        let exp = 3_u32.pow(s) * &t;

        // compute cubic root of (a^t)^-1, say h
        let mut h = ark_bn254::Fq12::ONE;
        let cc = c.pow([e as u64]);
        let mut c = c.inverse().unwrap();

        for i in 1..(s as i32) {
            let delta = (s as i32) - i - 1;
            let d = if delta < 0 {
                r.pow((&exp / 3_u32.pow((-delta) as u32)).to_u64_digits())
            } else {
                r.pow([3_u32.pow(delta as u32).to_u64().unwrap()])
            };
            if d == cc {
                (h, r) = (h * c, r * c.pow([3 as u64]));
            } else if d == cc.pow([2_u64]) {
                (h, r) = (h * c.pow([2 as u64]), r * c.pow([3 as u64]).pow([2 as u64]));
            }
            c = c.pow([3 as u64])
        }

        // recover cubic root of a
        // r = a^k * h
        r = a.pow(k.to_u64_digits()) * h;
        if t == 3_u32 * k + 1_u32 {
            r = r.inverse().unwrap();
        }

        assert_eq!(r.pow([3 as u64]), a);
        r
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::params::MODULUS;
    use std::ops::Deref;

    use ark_ff::{Field, One};
    use ark_std::UniformRand;
    use num_bigint::BigUint;

    use crate::params;

    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_compute_c_wi() {
        // 1. constant params

        let s = 3_u32;
        let exp = MODULUS.pow(12) - 1_u32;
        let h = &exp / params::R.deref();
        let t = &exp / 3_u32.pow(s);
        let k = (&t + 1_u32) / 3_u32;
        let m = params::LAMBDA.deref() / params::R.deref();
        let d = 3_u32;
        let mm = &m / d;

        let mut prng = ChaCha20Rng::seed_from_u64(0);
        let cofactor_cubic = 3_u32.pow(s - 1) * &t;

        // sample a miller loop result f which is cubic non-residue
        let f = {
            // (p^12 - 1) // 3
            let mut f = ark_bn254::Fq12::rand(&mut prng).pow(params::R.deref().to_u64_digits());
            let mut legendre = f.pow(cofactor_cubic.to_u64_digits());
            while legendre == ark_bn254::Fq12::ONE {
                f = ark_bn254::Fq12::rand(&mut prng).pow(params::R.deref().to_u64_digits());
                legendre = f.pow(cofactor_cubic.to_u64_digits());
            }
            f
        };
        assert_eq!(f.pow(h.to_u64_digits()), ark_bn254::Fq12::ONE);

        // sample a proper scalar w which is cubic non-residue
        let w = {
            let (mut w, mut z) = (ark_bn254::Fq12::one(), ark_bn254::Fq12::one());
            while w == ark_bn254::Fq12::ONE {
                // choose z which is 3-th non-residue
                let mut legendre = ark_bn254::Fq12::ONE;
                while legendre == ark_bn254::Fq12::ONE {
                    z = ark_bn254::Fq12::rand(&mut prng);
                    legendre = z.pow(cofactor_cubic.to_u64_digits());
                }
                // obtain w which is t-th power of z
                w = z.pow(t.to_u64_digits());
            }
            w
        };
        // make sure 27-th root w, is 3-th non-residue and r-th residue
        assert_ne!(w.pow(cofactor_cubic.to_u64_digits()), ark_bn254::Fq12::ONE);
        assert_eq!(w.pow(h.to_u64_digits()), ark_bn254::Fq12::ONE);

        // just two option, w and w^2, since w^3 must be cubic residue, leading f*w^3 must not be cubic residue
        let mut wi = w;
        if (f * wi).pow(cofactor_cubic.to_u64_digits()) != ark_bn254::Fq12::ONE {
            assert_eq!(
                (f * w * w).pow(cofactor_cubic.to_u64_digits()),
                ark_bn254::Fq12::ONE
            );
            wi = w * w;
        }
        assert_eq!(wi.pow(h.to_u64_digits()), ark_bn254::Fq12::ONE);

        assert_eq!(params::LAMBDA.clone(), &d * &mm * params::R.deref());
        // f1 is scaled f
        let f1 = f * wi;

        // r-th root of f1, say f2
        let r_inv = params::R.deref().modinv(&h).unwrap();
        assert_ne!(r_inv, BigUint::one());
        let f2 = f1.pow(r_inv.to_u64_digits());

        // m'-th root of f, say f3
        let mm_inv = mm.modinv(&(params::R.deref() * h)).unwrap();
        assert_ne!(mm_inv, BigUint::one());
        let f3 = f2.pow(mm_inv.to_u64_digits());

        // d-th (cubic) root, say c
        let c = LambdaResidues::tonelli_shanks_cubic(f3, w, s, t, k);
        assert_eq!(c.pow(params::LAMBDA.deref().to_u64_digits()), f * wi);
    }
}
