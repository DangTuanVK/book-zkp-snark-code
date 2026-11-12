//! Lab R1 (Chương 5): Groth16 (ark-groth16 0.6) và Marlin (ark-marlin 0.3) trên
//! BLS12-381 cho cùng một mạch, nhân rộng tới miền 2^10 .. 2^16.
//!
//! Mạch: lặp c bản sao của mạch ví dụ Chương 5. Bản sao thứ k có đầu vào công
//! khai N_k và witness a_k, b_k; mỗi bản sao gồm 13 ràng buộc R1CS:
//!   - 5 + 5 ràng buộc bit  (b_i * b_i = b_i),
//!   - 2 ràng buộc tổ hợp   (sum 2^i b_i) * 1 = a   (và tương tự cho b),
//!   - 1 ràng buộc nhân     a * b = N_k.
//! Chọn c = floor((2^k - 1) / 14) để Groth16 dùng miền đúng 2^k
//! (ark-groth16: cỡ miền = số ràng buộc + số biến instance, làm tròn lên 2^j).
//!
//! Cách chạy:  cargo run --release -- [k1,k2,...] [số_lần] [groth16|marlin|ca_hai]
//! Ví dụ:      RAYON_NUM_THREADS=1 cargo run --release -- 10,12 5 ca_hai
//! Mỗi dòng kết quả in dạng CSV (tiền tố "CSV,") để gom vào ket_qua/ch05/lab_r1.csv.

use std::env;
use std::time::{Duration, Instant};

/// Cặp (a_k, b_k) tất định, 2 <= a, b <= 31 (không tràn trong F_r của BLS12-381).
fn cap_thua_so(c: usize) -> Vec<(u64, u64)> {
    (0..c)
        .map(|k| (2 + (k as u64 % 30), 2 + ((7 * k as u64 + 3) % 30)))
        .collect()
}

fn so_ban_sao(k: u32) -> usize {
    ((1usize << k) - 1) / 14
}

fn trung_vi(mut v: Vec<Duration>) -> Duration {
    v.sort();
    v[v.len() / 2]
}

fn ms(d: Duration) -> f64 {
    d.as_secs_f64() * 1e3
}

fn tai_may() -> String {
    std::fs::read_to_string("/proc/loadavg")
        .map(|s| s.split_whitespace().take(3).collect::<Vec<_>>().join(" "))
        .unwrap_or_else(|_| "?".into())
}

// ======================= Groth16 trên hệ arkworks 0.6 =======================
mod g16 {
    use super::*;
    use ark_bls12_381::{Bls12_381, Fr};
    use ark_groth16::Groth16;
    use ark_relations::{
        gr1cs::{ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef, SynthesisError},
        lc,
    };
    use ark_serialize::CanonicalSerialize;
    use ark_snark::{CircuitSpecificSetupSNARK, SNARK};
    use ark_std::rand::{rngs::StdRng, SeedableRng};

    #[derive(Clone)]
    pub struct LapPhanTich {
        pub cap: Vec<(u64, u64)>,
        pub co_gia_tri: bool,
    }

    impl ConstraintSynthesizer<Fr> for LapPhanTich {
        fn generate_constraints(self, cs: ConstraintSystemRef<Fr>) -> Result<(), SynthesisError> {
            let cgt = self.co_gia_tri;
            let thieu = || SynthesisError::AssignmentMissing;
            for &(a, b) in &self.cap {
                let n = cs.new_input_variable(|| if cgt { Ok(Fr::from(a * b)) } else { Err(thieu()) })?;
                let mut bien = Vec::new();
                for &x in &[a, b] {
                    let v = cs.new_witness_variable(|| if cgt { Ok(Fr::from(x)) } else { Err(thieu()) })?;
                    let mut tong = lc!();
                    for i in 0..5 {
                        let bit = (x >> i) & 1;
                        let bi = cs.new_witness_variable(|| if cgt { Ok(Fr::from(bit)) } else { Err(thieu()) })?;
                        cs.enforce_r1cs_constraint(|| lc!() + bi, || lc!() + bi, || lc!() + bi)?;
                        tong = tong + (Fr::from(1u64 << i), bi);
                    }
                    let t = tong.clone();
                    cs.enforce_r1cs_constraint(move || t, || lc!() + ark_relations::gr1cs::Variable::One, || lc!() + v)?;
                    bien.push(v);
                }
                let (va, vb) = (bien[0], bien[1]);
                cs.enforce_r1cs_constraint(|| lc!() + va, || lc!() + vb, || lc!() + n)?;
            }
            Ok(())
        }
    }

    /// Trả về (số ràng buộc, số biến instance, số biến witness, thỏa?)
    pub fn thong_ke(m: &LapPhanTich) -> (usize, usize, usize, bool) {
        let cs = ConstraintSystem::<Fr>::new_ref();
        m.clone().generate_constraints(cs.clone()).unwrap();
        cs.finalize();
        (cs.num_constraints(), cs.num_instance_variables(), cs.num_witness_variables(),
         cs.is_satisfied().unwrap())
    }

    pub fn do_luong(k: u32, so_lan: usize) -> Vec<String> {
        let c = so_ban_sao(k);
        let cap = cap_thua_so(c);
        let mach = LapPhanTich { cap: cap.clone(), co_gia_tri: true };
        let rong = LapPhanTich { cap: cap.clone(), co_gia_tri: false };
        let (nc, ni, nw, ok) = thong_ke(&mach);
        assert!(ok);
        let mien = (nc + ni).next_power_of_two();
        let x: Vec<Fr> = cap.iter().map(|&(a, b)| Fr::from(a * b)).collect();
        let mut rng = StdRng::seed_from_u64(391 + k as u64);
        let tai_truoc = tai_may();

        let (mut t_setup, mut t_prove, mut t_verify) = (vec![], vec![], vec![]);
        let mut co_kich_thuoc = 0;
        let mut dung = true;
        for _ in 0..so_lan {
            let t = Instant::now();
            let (pk, vk) = Groth16::<Bls12_381>::setup(rong.clone(), &mut rng).unwrap();
            t_setup.push(t.elapsed());
            let pvk = Groth16::<Bls12_381>::process_vk(&vk).unwrap();
            let t = Instant::now();
            let proof = Groth16::<Bls12_381>::prove(&pk, mach.clone(), &mut rng).unwrap();
            t_prove.push(t.elapsed());
            let t = Instant::now();
            let v = Groth16::<Bls12_381>::verify_with_processed_vk(&pvk, &x, &proof).unwrap();
            t_verify.push(t.elapsed());
            dung &= v;
            co_kich_thuoc = proof.compressed_size();
            // Sửa một đầu vào công khai: phải bị bác
            let mut x2 = x.clone();
            x2[0] += Fr::from(1u64);
            dung &= !Groth16::<Bls12_381>::verify_with_processed_vk(&pvk, &x2, &proof).unwrap();
        }
        let tai_sau = tai_may();
        assert!(dung);
        vec![format!(
            "CSV,groth16,{k},{c},{nc},{ni},{nw},{mien},{so_lan},{:.3},{:.3},{:.3},{co_kich_thuoc},{dung},{tai_truoc},{tai_sau}",
            ms(trung_vi(t_setup)), ms(trung_vi(t_prove)), ms(trung_vi(t_verify))
        )]
    }
}

// ======================= Marlin trên hệ arkworks 0.3 ========================
mod marlin {
    use super::*;
    use ark_bls12_381_v03::{Bls12_381, Fr};
    use ark_marlin::Marlin;
    use ark_poly_commit_v03::marlin_pc::MarlinKZG10;
    use ark_poly_v03::univariate::DensePolynomial;
    use ark_relations_v03::{
        lc,
        r1cs::{ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef, OptimizationGoal,
               SynthesisError, SynthesisMode, Variable},
    };
    use ark_serialize_v03::CanonicalSerialize;
    use ark_std_v03::rand::{rngs::StdRng, SeedableRng};
    use blake2_v09::Blake2s;

    type Pc = MarlinKZG10<Bls12_381, DensePolynomial<Fr>>;
    type MarlinInst = Marlin<Fr, Pc, Blake2s>;

    #[derive(Clone)]
    pub struct LapPhanTich {
        pub cap: Vec<(u64, u64)>,
    }

    impl ConstraintSynthesizer<Fr> for LapPhanTich {
        fn generate_constraints(self, cs: ConstraintSystemRef<Fr>) -> Result<(), SynthesisError> {
            for &(a, b) in &self.cap {
                let n = cs.new_input_variable(|| Ok(Fr::from(a * b)))?;
                let mut bien = Vec::new();
                for &x in &[a, b] {
                    let v = cs.new_witness_variable(|| Ok(Fr::from(x)))?;
                    let mut tong = lc!();
                    for i in 0..5 {
                        let bit = (x >> i) & 1;
                        let bi = cs.new_witness_variable(|| Ok(Fr::from(bit)))?;
                        cs.enforce_constraint(lc!() + bi, lc!() + bi, lc!() + bi)?;
                        tong = tong + (Fr::from(1u64 << i), bi);
                    }
                    cs.enforce_constraint(tong, lc!() + Variable::One, lc!() + v)?;
                    bien.push(v);
                }
                cs.enforce_constraint(lc!() + bien[0], lc!() + bien[1], lc!() + n)?;
            }
            Ok(())
        }
    }

    /// Kích thước Marlin cần cho SRS: (cạnh ma trận sau khi đệm vuông, số phần tử khác 0)
    fn kich_thuoc_index(m: &LapPhanTich) -> (usize, usize, usize) {
        let cs = ConstraintSystem::<Fr>::new_ref();
        cs.set_optimization_goal(OptimizationGoal::Weight);
        cs.set_mode(SynthesisMode::Setup);
        m.clone().generate_constraints(cs.clone()).unwrap();
        cs.finalize();
        let mt = cs.to_matrices().unwrap();
        let nnz = mt.a_num_non_zero.max(mt.b_num_non_zero).max(mt.c_num_non_zero);
        let ni = cs.num_instance_variables().next_power_of_two();
        let nv = ni + cs.num_witness_variables();
        (nv.max(cs.num_constraints()), nnz, cs.num_constraints())
    }

    pub fn do_luong(k: u32, so_lan: usize) -> Vec<String> {
        let c = so_ban_sao(k);
        let cap = cap_thua_so(c);
        let mach = LapPhanTich { cap: cap.clone() };
        let (canh, nnz, nc) = kich_thuoc_index(&mach);
        let x: Vec<Fr> = cap.iter().map(|&(a, b)| Fr::from(a * b)).collect();
        let mut rng = StdRng::seed_from_u64(449 + k as u64);
        let tai_truoc = tai_may();

        let t = Instant::now();
        let srs = MarlinInst::universal_setup(canh, canh, nnz, &mut rng).unwrap();
        let t_srs = t.elapsed();

        let (mut t_index, mut t_prove, mut t_verify) = (vec![], vec![], vec![]);
        let mut co = 0;
        let mut dung = true;
        for _ in 0..so_lan {
            let t = Instant::now();
            let (ipk, ivk) = MarlinInst::index(&srs, mach.clone()).unwrap();
            t_index.push(t.elapsed());
            let t = Instant::now();
            let proof = MarlinInst::prove(&ipk, mach.clone(), &mut rng).unwrap();
            t_prove.push(t.elapsed());
            let t = Instant::now();
            let v = MarlinInst::verify(&ivk, &x, &proof, &mut rng).unwrap();
            t_verify.push(t.elapsed());
            dung &= v;
            co = proof.serialized_size();
            let mut x2 = x.clone();
            x2[0] += Fr::from(1u64);
            dung &= !MarlinInst::verify(&ivk, &x2, &proof, &mut rng).unwrap();
        }
        let tai_sau = tai_may();
        assert!(dung);
        vec![format!(
            "CSV,marlin,{k},{c},{nc},{},{},{},{so_lan},{:.3},{:.3},{:.3},{co},{dung},{tai_truoc},{tai_sau},srs_ms={:.3},canh={canh},nnz={nnz},mien_K={}",
            c + 1, 12 * c, canh.next_power_of_two(),
            ms(trung_vi(t_index)), ms(trung_vi(t_prove)), ms(trung_vi(t_verify)), ms(t_srs), nnz.next_power_of_two()
        )]
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let ks: Vec<u32> = args.get(1).map(|s| s.split(',').map(|x| x.parse().unwrap()).collect())
        .unwrap_or(vec![10, 12, 14, 16]);
    let so_lan: usize = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(5);
    let he = args.get(3).cloned().unwrap_or("ca_hai".into());
    let luong = env::var("RAYON_NUM_THREADS").unwrap_or("mặc định (= số CPU logic)".into());
    println!("# Lab R1: RAYON_NUM_THREADS = {luong}; số lần đo mỗi điểm = {so_lan} (lấy trung vị)");
    println!("# cột: he,k,ban_sao,rang_buoc,bien_instance,bien_witness,mien,so_lan,setup_ms,prove_ms,verify_ms,proof_byte,dung,tai_truoc,tai_sau[,phu]");
    println!("# Groth16: setup = sinh CRS riêng cho mạch. Marlin: setup = index (chuyên biệt hóa từ SRS phổ quát);");
    println!("#          srs_ms = thời gian universal_setup (một lần); rang_buoc Marlin = trước khi đệm vuông; mien = |H| (cạnh ma trận sau đệm, làm tròn 2^j), mien_K = |K|");
    for &k in &ks {
        if he == "groth16" || he == "ca_hai" {
            for l in g16::do_luong(k, so_lan) { println!("{l}"); }
        }
        if he == "marlin" || he == "ca_hai" {
            for l in marlin::do_luong(k, so_lan) { println!("{l}"); }
        }
    }
}
