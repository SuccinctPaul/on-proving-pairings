use crate::miller_lines::MillerLines;
use crate::params;
use crate::params::MODULUS;
use crate::utils::biguint_to_naf;
use ark_bn254::{Fq, Fq12, Fq2, Fq6, G1Projective, G2Projective};
use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::Field;
use std::ops::{Add, Mul, Neg, Sub};

pub struct NativeMillerLoop;

impl NativeMillerLoop {
    fn line_func_add(
        r: G2Projective,
        p: G2Projective,
        q: G1Projective,
        r2: Fq2,
    ) -> (Fq2, Fq2, Fq2, G2Projective) {
        let r_t = r.z.square();
        let B = p.x * r_t;
        let mut D = p.y + r.z;
        D = D.square();
        D -= r2;
        D -= r_t;
        D *= r_t;

        let H = B - r.x;
        let I = H.square();

        let E = I.double().double();

        let J = H * E;
        let mut L1 = D - r.y;
        L1 -= r.y;

        let V = r.x * E;

        let mut r_x = L1.square();
        r_x -= J;
        r_x -= V.double();

        let mut r_z = r.z + H;
        r_z = r_z.square();
        r_z -= r_t;
        r_z -= I;

        let mut t = V - r_x;
        t *= L1;
        let mut t2 = r.y * J;
        t2 = t2.double();
        let r_y = t - t2;

        let r_out = G2Projective::new(r_x, r_y, r_z);

        let mut t = p.y + r_z;
        t = t.square();
        t = t - r2;
        t = t - (r_z.square());

        let mut t2 = L1 * p.x;
        t2 = t2.double();

        let a = t2 - t;

        let x = q.y;
        let mut c = r_z;
        c.mul_assign_by_basefield(&q.y);
        c = c.double();

        let mut b = L1.neg();
        b.mul_assign_by_basefield(&q.x);
        b = b.double();

        // abandon the convenience of projective coordinate, be consistent with verifier
        // 2 * z_r
        let aux_inv = r_z.double().inverse().unwrap();

        let a = a.mul(aux_inv);
        let b = b.mul(aux_inv);
        let c = c.mul(aux_inv);

        (a, b, c, r_out)
    }

    fn line_func_double(r: G2Projective, q: G1Projective) -> (Fq2, Fq2, Fq2, G2Projective) {
        let r_t = r.z.square();

        let A = r.x.square();
        let B = r.y.square();
        let C = B.square();

        let mut D = r.x + B;
        D = D.square();
        D -= A;
        D -= C;
        D = D.double();

        let E = A.double() + A;
        let F = E.square();

        // C^8
        let C8 = C.double().double().double();

        let r_x = F - D.double();
        let r_y = E * (D - r_x) - C8;

        // (y+z)*(y+z) - (y*y) - (z*z) = 2*y*z
        let r_z = (r.y + r.z).square() - B - r_t;

        assert_eq!(r_z, r.y * r.z.double());

        let r_out = G2Projective::new(r_x, r_y, r_z);
        assert!(r_out.into_affine().is_on_curve());

        let mut a = r.x + E;
        a = a.square();
        a -= A + F + B.double().double();

        let mut t = E * r_t;
        t = t.double();
        let mut b = t.neg();
        b.mul_assign_by_basefield(&q.x);

        let mut c = r_z * r_t;
        c = c.double();
        c.mul_assign_by_basefield(&q.y);

        // abandon the convenience of projective coordinate, be consistant with verifier
        // 2 * z_r * z_t^2
        let aux_inv = r_t.mul(r_z).double().inverse().unwrap();
        let a = a.mul(aux_inv);
        let b = b.mul(aux_inv);
        let c = c.mul(aux_inv);

        (a, b, c, r_out)
    }

    pub fn mul_line(r: Fq12, a: Fq2, b: Fq2, c: Fq2) -> Fq12 {
        // See function fp12e_mul_line in dclxvi
        let t1 = Fq6::new(Fq2::ZERO, a, b);
        let t2 = Fq6::new(Fq2::ZERO, a, b.add(c));

        let t1 = t1.mul(r.c0);
        let mut t3 = r.c1;
        t3.mul_by_fp2(&c);

        let y = t3;
        let mut x = r.c0.add(r.c1);
        x = x.mul(t2);
        x = x.sub(t1);
        let x = x.sub(y);
        // TODO: y = y.add(t1.mul_tau())
        // let y = y.add()

        Fq12::new(x, y)
    }

    // Native miller loop: With addend of p^3
    pub fn miller_loop(p: G1Projective, q: G2Projective) -> Fq12 {
        let P = p.into_affine();
        let Q = q.into_affine();
        // let P = p;
        // let Q = q;
        let mQ = Q.neg();

        // 6x + 2 in NAF
        let mut naf_digits = biguint_to_naf(params::E.clone());
        naf_digits.reverse();
        naf_digits.remove(0);

        let mut T = q;
        let Qp = Q.y.square();

        let mut f = Fq12::new(Fq6::ONE, Fq6::ZERO);
        // f_list is used for debug
        let mut f_list = vec![];

        naf_digits.iter().enumerate().for_each(|(i, digit)| {
            f = f.square();
            let (a, b, c, t) = Self::line_func_double(G2Projective::from(T), G1Projective::from(P));
            T = t;
            f = MillerLines::mul_line_base(f, a, b, c);
            f_list.push(f);

            if *digit == 1 {
                let (a, b, c, t) = Self::line_func_add(
                    G2Projective::from(T),
                    G2Projective::from(Q),
                    G1Projective::from(P),
                    Qp,
                );
                T = t;

                f = MillerLines::mul_line_base(f, a, b, c);
                f_list.push(f);
            } else if *digit == -1 {
                let (a, b, c, t) = Self::line_func_add(
                    G2Projective::from(T),
                    G2Projective::from(mQ),
                    G1Projective::from(P),
                    Qp,
                );
                T = t;

                f = MillerLines::mul_line_base(f, a, b, c);

                f_list.push(f);
            }
        });

        assert_eq!(T, Q.mul_bigint(params::E.to_u64_digits()));

        // aaaa
        let (mut x, mut y) = (Q.x.clone(), Q.y.clone());

        let pi_1_Q = G2Projective::new(
            x.conjugate_in_place().mul(params::BETA_PI_1[1]),
            y.conjugate_in_place().mul(params::BETA_PI_1[2]),
            Fq2::ONE,
        );
        assert_eq!(pi_1_Q, Q.mul_bigint(params::MODULUS.to_u64_digits()));

        // 2.2. Q2 = pi2(Q)
        // x = x * beta * (2 * (p^2 - 1) / 6)
        // y = y * beta * (3 * (p^2 - 1) / 6) = -y
        let (x, y) = (Q.x, Q.y);
        let pi_2_Q = G2Projective::new(
            x.mul(params::BETA_PI_2[1]),
            y.mul(params::BETA_PI_2[2]),
            Fq2::ONE,
        );
        assert_eq!(pi_2_Q, Q.mul_bigint(MODULUS.pow(2).to_u64_digits()));

        // 2.3. Q3 = pi3(Q)
        // x = x' * beta * (2 * (p^3 - 1) / 6)
        // y = y' * beta * (3 * (p^3 - 1) / 6)
        let (mut x, mut y) = (Q.x.clone(), Q.y.clone());

        let pi_3_Q = G2Projective::new(
            x.conjugate_in_place().mul(params::BETA_PI_3[1]),
            y.conjugate_in_place().mul(params::BETA_PI_3[2]),
            Fq2::ONE,
        );
        assert_eq!(pi_3_Q, Q.mul_bigint(MODULUS.pow(3).to_u64_digits()));

        let Qp = pi_1_Q.y.square();
        let (a, b, c, t) = Self::line_func_add(
            G2Projective::from(T),
            G2Projective::from(pi_1_Q),
            G1Projective::from(P),
            Qp,
        );
        T = t;

        let f = MillerLines::mul_line_base(f, a, b, c);
        f_list.push(f);

        let Qp = pi_2_Q.y.square();
        // let (a, b, c, T) = line_func_add(T, pi_2_Q.neg(), P, Qp);
        let (a, b, c, t) = Self::line_func_add(
            G2Projective::from(T),
            G2Projective::from(pi_2_Q.neg()),
            G1Projective::from(P),
            Qp,
        );
        T = t;

        let f = MillerLines::mul_line_base(f, a, b, c);
        f_list.push(f);

        // k = 6 * x + 2 + px(x) - px(x) ** 2
        // assert(T == Q.scalar_mul(k if k > 0 else rx(x) - ((-k) % rx(x))))
        // assert(T.is_infinite() == False)

        // The bn254's miller donest's have the following
        let eval = Fq12::new(
            Fq6::new(
                Fq2::new(P.x, Fq::ZERO),
                T.x.mul(T.z.inverse().unwrap().square()).neg(),
                Fq2::ZERO,
            ),
            Fq6::ZERO,
        );

        let T = T.add(pi_3_Q);
        let f = f.mul(eval);
        f_list.push(f);
        // k = 6 * x + 2 + px(x) - px(x) ** 2 + px(x) ** 3
        assert_eq!(T, Q.mul_bigint(params::LAMBDA.to_u64_digits()));
        // assert(T.is_infinite() == True)

        f
    }
}

#[cfg(test)]
mod test {
    use super::*;
<<<<<<<< HEAD:on-proving-pairings/src/optimal_ate.rs
    use crate::dev;
========
    use crate::constants;
>>>>>>>> d389fd5 (feat:  add setup/prove/verify concept and refactor the project structure  (#5)):poc/src/optimal_ate.rs
    use ark_bn254::{Bn254, G1Affine, G2Affine};
    use ark_ec::pairing::Pairing;
    use ark_ec::{AffineRepr, CurveGroup};
    use num_bigint::BigUint;
    use num_traits::FromPrimitive;

    #[test]
    fn test_bn254_miller_loop() {
        let Q = dev::g2
            .mul_bigint(BigUint::from_i8(3).unwrap().to_u64_digits())
            .into_affine();
        let P = dev::g1
            .mul_bigint(BigUint::from_i8(4).unwrap().to_u64_digits())
            .into_affine();

        println!("P:{:?}", P);
        println!("\n Q:{:?}", Q);
        let actual = Bn254::miller_loop(
            <G1Affine as Into<<Bn254 as Pairing>::G1Prepared>>::into(P),
            <G2Affine as Into<<Bn254 as Pairing>::G2Prepared>>::into(Q),
        );

        println!("\n actual: {:?}", actual.0.to_string());
    }

    #[test]
    fn test_native_miller_loop() {
        let Q = dev::g2.mul_bigint(BigUint::from_i8(3).unwrap().to_u64_digits());
        let P = dev::g1.mul_bigint(BigUint::from_i8(4).unwrap().to_u64_digits());

        println!("P:{:?}", P);
        println!("\n Q:{:?}", Q);
        let actual = NativeMillerLoop::miller_loop(P, Q);

        println!("\n actual: {:?}", actual.to_string());
    }
}
