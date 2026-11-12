//! Lab R1 (bản chạy thử): Groth16 trên BLS12-381 cho mạch R1CS tí hon a·b = c.
//! Phần 2: Marlin (SNARK phổ quát, SRS dùng chung) qua ark-marlin 0.3 (hệ ark 0.3).
//! Bài lab đầy đủ viết ở đợt sau.

use ark_bls12_381::{Bls12_381, Fr};
use ark_groth16::Groth16;
use ark_relations::{
    gr1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError},
    lc,
};
use ark_serialize::CanonicalSerialize;
use ark_snark::{CircuitSpecificSetupSNARK, SNARK};
use ark_std::rand::{rngs::StdRng, SeedableRng};
use std::time::Instant;

// Mạch: biết a, b (bí mật) sao cho a·b = c (công khai)
struct NhanHaiSo {
    a: Option<Fr>,
    b: Option<Fr>,
}

impl ConstraintSynthesizer<Fr> for NhanHaiSo {
    fn generate_constraints(self, cs: ConstraintSystemRef<Fr>) -> Result<(), SynthesisError> {
        let a = cs.new_witness_variable(|| self.a.ok_or(SynthesisError::AssignmentMissing))?;
        let b = cs.new_witness_variable(|| self.b.ok_or(SynthesisError::AssignmentMissing))?;
        let c = cs.new_input_variable(|| {
            let a = self.a.ok_or(SynthesisError::AssignmentMissing)?;
            let b = self.b.ok_or(SynthesisError::AssignmentMissing)?;
            Ok(a * b)
        })?;
        // Một ràng buộc R1CS duy nhất: a * b = c
        cs.enforce_r1cs_constraint(|| lc!() + a, || lc!() + b, || lc!() + c)?;
        Ok(())
    }
}

fn main() {
    // Seed cố định để chạy lại cho cùng kết quả
    let mut rng = StdRng::seed_from_u64(391);

    let t = Instant::now();
    let (pk, vk) = Groth16::<Bls12_381>::setup(NhanHaiSo { a: None, b: None }, &mut rng)
        .expect("setup lỗi");
    println!("[R1] Groth16/BLS12-381 setup: {:.2?}", t.elapsed());

    let (a, b) = (Fr::from(17u64), Fr::from(23u64));
    let c = Fr::from(391u64);

    let t = Instant::now();
    let proof = Groth16::<Bls12_381>::prove(&pk, NhanHaiSo { a: Some(a), b: Some(b) }, &mut rng)
        .expect("prove lỗi");
    println!("[R1] prove (a=17, b=23): {:.2?}", t.elapsed());
    println!("[R1] kích thước proof (nén): {} byte", proof.compressed_size());

    let t = Instant::now();
    let ok = Groth16::<Bls12_381>::verify(&vk, &[c], &proof).expect("verify lỗi");
    println!("[R1] verify với c=391: {} ({:.2?})", ok, t.elapsed());

    // Đầu vào công khai sai phải bị bác
    let sai = Groth16::<Bls12_381>::verify(&vk, &[Fr::from(392u64)], &proof).unwrap();
    println!("[R1] verify với c=392: {}", sai);
    assert!(ok && !sai);

    marlin::chay();
}

// ---- Marlin trên hệ arkworks 0.3 (kiểu dữ liệu không trộn được với ark 0.6) ----
mod marlin {
    use ark_bls12_381_v03::{Bls12_381, Fr};
    use ark_marlin::Marlin;
    use ark_poly_commit_v03::marlin_pc::MarlinKZG10;
    use ark_poly_v03::univariate::DensePolynomial;
    use ark_relations_v03::{
        lc,
        r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError},
    };
    use blake2_v09::Blake2s;
    use std::time::Instant;

    type Pc = MarlinKZG10<Bls12_381, DensePolynomial<Fr>>;
    type MarlinInst = Marlin<Fr, Pc, Blake2s>;
    const SO_LAN_LAP: usize = 2; // = 1 thì panic trong ark-poly-commit 0.3

    // Cùng mạch a·b = c, viết lại bằng API ark-relations 0.3
    #[derive(Clone, Copy)]
    struct NhanHaiSo {
        a: Option<Fr>,
        b: Option<Fr>,
    }

    impl ConstraintSynthesizer<Fr> for NhanHaiSo {
        fn generate_constraints(self, cs: ConstraintSystemRef<Fr>) -> Result<(), SynthesisError> {
            let a = cs.new_witness_variable(|| self.a.ok_or(SynthesisError::AssignmentMissing))?;
            let b = cs.new_witness_variable(|| self.b.ok_or(SynthesisError::AssignmentMissing))?;
            let c = cs.new_input_variable(|| {
                let a = self.a.ok_or(SynthesisError::AssignmentMissing)?;
                let b = self.b.ok_or(SynthesisError::AssignmentMissing)?;
                Ok(a * b)
            })?;
            // ark-marlin 0.3 lỗi (panic) với mạch quá nhỏ: lặp ràng buộc để
            // miền H, K đủ lớn; ý nghĩa mạch không đổi
            for _ in 0..SO_LAN_LAP {
                cs.enforce_constraint(lc!() + a, lc!() + b, lc!() + c)?;
            }
            Ok(())
        }
    }

    pub fn chay() {
        let rng = &mut ark_std_v03::test_rng();

        // SRS phổ quát: tối đa 64 ràng buộc, 64 biến, 64 phần tử khác 0
        // (SRS quá sát cỡ mạch gây panic trong ark-poly-commit 0.3)
        let t = Instant::now();
        let srs = MarlinInst::universal_setup(64, 64, 64, rng).expect("universal_setup lỗi");
        println!("[R1] Marlin/BLS12-381 universal setup (64,64,64): {:.2?}", t.elapsed());

        let mach = NhanHaiSo { a: Some(Fr::from(17u64)), b: Some(Fr::from(23u64)) };
        let t = Instant::now();
        let (ipk, ivk) = MarlinInst::index(&srs, mach).expect("index lỗi");
        println!("[R1] Marlin index (chuyên biệt hóa từ SRS): {:.2?}", t.elapsed());

        let t = Instant::now();
        let proof = MarlinInst::prove(&ipk, mach, rng).expect("prove lỗi");
        println!("[R1] Marlin prove (a=17, b=23): {:.2?}", t.elapsed());

        let t = Instant::now();
        let ok = MarlinInst::verify(&ivk, &[Fr::from(391u64)], &proof, rng).expect("verify lỗi");
        println!("[R1] Marlin verify với c=391: {} ({:.2?})", ok, t.elapsed());
        let sai = MarlinInst::verify(&ivk, &[Fr::from(392u64)], &proof, rng).expect("verify lỗi");
        println!("[R1] Marlin verify với c=392: {}", sai);
        assert!(ok && !sai);
    }
}
