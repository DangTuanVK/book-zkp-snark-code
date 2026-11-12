//! Lab R2 (bản chạy thử): IVC bằng Nova (nova-snark) trên chu trình Pallas/Vesta.
//! Mạch bước F: z_{i+1} = z_i^2 + 3 (một ràng buộc R1CS / bước).
//! Sau IVC, nén proof bằng Spartan + IPA (không cần trusted setup).
//! Bài lab đầy đủ viết ở đợt sau.

use ff::PrimeField;
use nova_snark::{
    frontend::{num::AllocatedNum, ConstraintSystem, SynthesisError},
    nova::{CompressedSNARK, PublicParams, RecursiveSNARK},
    provider::{PallasEngine, VestaEngine},
    traits::{circuit::StepCircuit, snark::RelaxedR1CSSNARKTrait, Engine},
};
use std::time::Instant;

type E1 = PallasEngine;
type E2 = VestaEngine;
type F1 = <E1 as Engine>::Scalar;
type EE1 = nova_snark::provider::ipa_pc::EvaluationEngine<E1>;
type EE2 = nova_snark::provider::ipa_pc::EvaluationEngine<E2>;
type S1 = nova_snark::spartan::snark::RelaxedR1CSSNARK<E1, EE1>;
type S2 = nova_snark::spartan::snark::RelaxedR1CSSNARK<E2, EE2>;

const SO_BUOC: usize = 5;

// Mạch bước: nhận z_i, trả z_{i+1} = z_i^2 + 3
#[derive(Clone, Debug, Default)]
struct BinhPhuongCong3;

impl<F: PrimeField> StepCircuit<F> for BinhPhuongCong3 {
    fn arity(&self) -> usize {
        1
    }

    fn synthesize<CS: ConstraintSystem<F>>(
        &self,
        cs: &mut CS,
        z: &[AllocatedNum<F>],
    ) -> Result<Vec<AllocatedNum<F>>, SynthesisError> {
        let x = &z[0];
        // y = x^2 + 3, ràng buộc: x * x = y - 3
        let y = AllocatedNum::alloc(cs.namespace(|| "y"), || {
            let v = x.get_value().ok_or(SynthesisError::AssignmentMissing)?;
            Ok(v * v + F::from(3u64))
        })?;
        cs.enforce(
            || "x*x = y - 3",
            |lc| lc + x.get_variable(),
            |lc| lc + x.get_variable(),
            |lc| lc + y.get_variable() - (F::from(3u64), CS::one()),
        );
        Ok(vec![y])
    }
}

fn main() {
    let mach = BinhPhuongCong3;

    let t = Instant::now();
    let pp = PublicParams::<E1, E2, BinhPhuongCong3>::setup(&mach, &*S1::ck_floor(), &*S2::ck_floor())
        .expect("setup lỗi");
    println!("[R2] Nova PublicParams::setup (Pallas/Vesta): {:.2?}", t.elapsed());
    println!(
        "[R2] số ràng buộc mỗi bước (mạch chính, mạch phụ): {:?}",
        pp.num_constraints()
    );

    let z0 = vec![F1::from(2u64)];
    let t = Instant::now();
    let mut rs = RecursiveSNARK::<E1, E2, BinhPhuongCong3>::new(&pp, &mach, &z0).expect("new lỗi");
    for i in 0..SO_BUOC {
        let ti = Instant::now();
        rs.prove_step(&pp, &mach).expect("prove_step lỗi");
        println!("[R2] prove_step {}: {:.2?}", i, ti.elapsed());
    }
    println!("[R2] tổng thời gian IVC {} bước: {:.2?}", SO_BUOC, t.elapsed());

    let t = Instant::now();
    let zn = rs.verify(&pp, SO_BUOC, &z0).expect("RecursiveSNARK::verify lỗi");
    println!("[R2] RecursiveSNARK::verify: OK ({:.2?})", t.elapsed());

    // Đối chiếu với tính trực tiếp ngoài mạch
    let mut z = F1::from(2u64);
    for _ in 0..SO_BUOC {
        z = z.square() + F1::from(3u64);
    }
    println!("[R2] z_{} khớp tính trực tiếp: {}", SO_BUOC, zn[0] == z);
    // In 8 byte thấp của z_5 (biểu diễn little-endian) để đối chiếu
    let repr = zn[0].to_repr();
    let thap = u64::from_le_bytes(repr.as_ref()[..8].try_into().unwrap());
    println!("[R2] z_{} (64 bit thấp) = {}", SO_BUOC, thap);

    // Nén proof IVC bằng Spartan (RelaxedR1CSSNARK) + IPA
    let t = Instant::now();
    let (pk, vk) = CompressedSNARK::<_, _, _, S1, S2>::setup(&pp).expect("setup nén lỗi");
    println!("[R2] CompressedSNARK::setup: {:.2?}", t.elapsed());
    let t = Instant::now();
    let cs = CompressedSNARK::<_, _, _, S1, S2>::prove(&pp, &pk, &rs).expect("prove nén lỗi");
    println!("[R2] CompressedSNARK::prove: {:.2?}", t.elapsed());
    let bytes = bincode::serde::encode_to_vec(&cs, bincode::config::legacy()).expect("serialize lỗi");
    println!("[R2] kích thước CompressedSNARK (bincode): {} byte", bytes.len());
    let t = Instant::now();
    let ok = cs.verify(&vk, SO_BUOC, &z0).is_ok();
    println!("[R2] CompressedSNARK::verify: {} ({:.2?})", ok, t.elapsed());
    assert!(ok && zn[0] == z);
}
