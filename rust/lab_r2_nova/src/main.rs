//! Lab R2 (Chương 9): IVC bằng Nova (nova-snark 0.76) trên chu trình Pallas/Vesta.
//!
//! Hàm bước F ("chuỗi kiểm phân tích"): trạng thái z = (N, k).
//!   - witness bước i: (a_i, b_i), luân phiên (17, 23) và (23, 17);
//!   - ràng buộc: a_i, b_i vừa 5 bit (5 ràng buộc bit + 1 tổ hợp, mỗi số),
//!     a_i · b_i = N (1 ràng buộc) -- đúng 13 ràng buộc như mạch R1CS của sách;
//!   - đầu ra z' = (N, k + 1) (1 ràng buộc cho bộ đếm).
//! Tùy chọn: thêm `them` ràng buộc bình phương giả (u <- u^2) để đo thời gian
//! mỗi bước theo cỡ mạch bước.
//! Nén proof IVC bằng Spartan (RelaxedR1CSSNARK) + IPA, không cần trusted setup.
//!
//! Cách chạy: RAYON_NUM_THREADS=16 cargo run --release [-- nhanh]

use ff::{Field, PrimeField};
use nova_snark::{
    frontend::{num::AllocatedNum, ConstraintSystem, SynthesisError},
    nova::{CompressedSNARK, PublicParams, RecursiveSNARK},
    provider::{PallasEngine, VestaEngine},
    traits::{circuit::StepCircuit, snark::RelaxedR1CSSNARKTrait, Engine},
};
use std::io::Write;
use std::time::{Duration, Instant};

type E1 = PallasEngine;
type E2 = VestaEngine;
type F1 = <E1 as Engine>::Scalar;
type EE1 = nova_snark::provider::ipa_pc::EvaluationEngine<E1>;
type EE2 = nova_snark::provider::ipa_pc::EvaluationEngine<E2>;
type S1 = nova_snark::spartan::snark::RelaxedR1CSSNARK<E1, EE1>;
type S2 = nova_snark::spartan::snark::RelaxedR1CSSNARK<E2, EE2>;

const N_VD: u64 = 391;
const SO_BIT: usize = 5;

#[derive(Clone, Debug)]
struct KiemPhanTich {
    a: u64,
    b: u64,
    rong: bool, // true: mạch bước rỗng (chỉ trả lại z), để đo phần mạch xác minh của Nova
    them: usize, // số ràng buộc bình phương thêm
}

impl KiemPhanTich {
    fn buoc(i: usize, them: usize) -> Self {
        let (a, b) = if i % 2 == 0 { (17, 23) } else { (23, 17) };
        KiemPhanTich { a, b, rong: false, them }
    }
}

fn phan_tich_bit<F: PrimeField, CS: ConstraintSystem<F>>(
    cs: &mut CS,
    ten: &str,
    v: u64,
) -> Result<AllocatedNum<F>, SynthesisError> {
    let x = AllocatedNum::alloc(cs.namespace(|| ten.to_string()), || Ok(F::from(v)))?;
    let mut bits = Vec::new();
    for j in 0..SO_BIT {
        let bj = AllocatedNum::alloc(cs.namespace(|| format!("{ten}_{j}")), || {
            Ok(F::from((v >> j) & 1))
        })?;
        // b·b = b
        cs.enforce(
            || format!("bit {ten}_{j}"),
            |lc| lc + bj.get_variable(),
            |lc| lc + bj.get_variable(),
            |lc| lc + bj.get_variable(),
        );
        bits.push(bj);
    }
    // (Σ 2^j b_j)·1 = x
    cs.enforce(
        || format!("tong {ten}"),
        |mut lc| {
            for (j, bj) in bits.iter().enumerate() {
                lc = lc + (F::from(1u64 << j), bj.get_variable());
            }
            lc
        },
        |lc| lc + CS::one(),
        |lc| lc + x.get_variable(),
    );
    Ok(x)
}

impl<F: PrimeField> StepCircuit<F> for KiemPhanTich {
    fn arity(&self) -> usize {
        2
    }

    fn synthesize<CS: ConstraintSystem<F>>(
        &self,
        cs: &mut CS,
        z: &[AllocatedNum<F>],
    ) -> Result<Vec<AllocatedNum<F>>, SynthesisError> {
        if self.rong {
            return Ok(vec![z[0].clone(), z[1].clone()]);
        }
        let (n, k) = (&z[0], &z[1]);
        let a = phan_tich_bit(cs, "a", self.a)?;
        let b = phan_tich_bit(cs, "b", self.b)?;
        // a·b = N
        cs.enforce(
            || "a*b = N",
            |lc| lc + a.get_variable(),
            |lc| lc + b.get_variable(),
            |lc| lc + n.get_variable(),
        );
        // bộ đếm k' = k + 1
        let k1 = AllocatedNum::alloc(cs.namespace(|| "k1"), || {
            Ok(k.get_value().ok_or(SynthesisError::AssignmentMissing)? + F::ONE)
        })?;
        cs.enforce(
            || "k' = k + 1",
            |lc| lc + k.get_variable() + CS::one(),
            |lc| lc + CS::one(),
            |lc| lc + k1.get_variable(),
        );
        // ràng buộc thêm: u_{j+1} = u_j^2, u_0 = a
        let mut u = a.clone();
        for j in 0..self.them {
            u = u.square(cs.namespace(|| format!("sq{j}")))?;
        }
        Ok(vec![n.clone(), k1])
    }
}

fn tai_may() -> String {
    std::fs::read_to_string("/proc/loadavg")
        .map(|s| s.split_whitespace().take(3).collect::<Vec<_>>().join(" "))
        .unwrap_or_default()
}

fn ms(d: Duration) -> f64 {
    d.as_secs_f64() * 1e3
}

fn trung_vi(v: &mut Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 1 { v[n / 2] } else { (v[n / 2 - 1] + v[n / 2]) / 2.0 }
}

struct KetQua {
    them: usize,
    so_buoc: usize,
    cons: (usize, usize),
    vars: (usize, usize),
    setup_ms: f64,
    buoc_tb: f64,
    buoc_tv: f64,
    buoc_max: f64,
    tong_ivc_ms: f64,
    verify_ms: f64,
    nen_setup_ms: f64,
    nen_prove_ms: f64,
    nen_verify_ms: f64,
    nen_byte: usize,
    rs_byte: usize,
    tai_truoc: String,
    tai_sau: String,
}

fn chay(so_buoc: usize, them: usize, nen: bool) -> KetQua {
    let tai_truoc = tai_may();
    let mau = KiemPhanTich::buoc(0, them);
    let t = Instant::now();
    let pp = PublicParams::<E1, E2, KiemPhanTich>::setup(&mau, &*S1::ck_floor(), &*S2::ck_floor())
        .expect("setup lỗi");
    let setup_ms = ms(t.elapsed());

    let z0 = vec![F1::from(N_VD), F1::ZERO];
    let t = Instant::now();
    let mut rs = RecursiveSNARK::<E1, E2, KiemPhanTich>::new(&pp, &KiemPhanTich::buoc(0, them), &z0)
        .expect("new lỗi");
    let mut tg = Vec::new();
    for i in 0..so_buoc {
        let c = KiemPhanTich::buoc(i, them);
        let ti = Instant::now();
        rs.prove_step(&pp, &c).expect("prove_step lỗi");
        if i > 0 {
            tg.push(ms(ti.elapsed())); // bước 0 chỉ khởi tạo (đã làm trong new)
        }
    }
    let tong_ivc_ms = ms(t.elapsed());
    let t = Instant::now();
    let zn = rs.verify(&pp, so_buoc, &z0).expect("RecursiveSNARK::verify lỗi");
    let verify_ms = ms(t.elapsed());
    assert_eq!(zn[0], F1::from(N_VD));
    assert_eq!(zn[1], F1::from(so_buoc as u64));
    let rs_byte = bincode::serde::encode_to_vec(&rs, bincode::config::legacy()).unwrap().len();

    let (mut nen_setup_ms, mut nen_prove_ms, mut nen_verify_ms, mut nen_byte) = (0.0, 0.0, 0.0, 0);
    if nen {
        let t = Instant::now();
        let (pk, vk) = CompressedSNARK::<_, _, _, S1, S2>::setup(&pp).expect("setup nén lỗi");
        nen_setup_ms = ms(t.elapsed());
        let t = Instant::now();
        let cs = CompressedSNARK::<_, _, _, S1, S2>::prove(&pp, &pk, &rs).expect("prove nén lỗi");
        nen_prove_ms = ms(t.elapsed());
        nen_byte = bincode::serde::encode_to_vec(&cs, bincode::config::legacy()).unwrap().len();
        let t = Instant::now();
        let ok = cs.verify(&vk, so_buoc, &z0);
        nen_verify_ms = ms(t.elapsed());
        assert!(ok.is_ok(), "CompressedSNARK::verify thất bại");
        // sai số bước -> phải từ chối
        assert!(cs.verify(&vk, so_buoc + 1, &z0).is_err());
    }
    let buoc_tb = tg.iter().sum::<f64>() / tg.len() as f64;
    let buoc_max = tg.iter().cloned().fold(0.0, f64::max);
    let buoc_tv = trung_vi(&mut tg);
    KetQua {
        them, so_buoc, cons: pp.num_constraints(), vars: pp.num_variables(), setup_ms,
        buoc_tb, buoc_tv, buoc_max, tong_ivc_ms, verify_ms, nen_setup_ms, nen_prove_ms,
        nen_verify_ms, nen_byte, rs_byte, tai_truoc, tai_sau: tai_may(),
    }
}

fn main() {
    let nhanh = std::env::args().any(|a| a == "nhanh");
    let luong = std::env::var("RAYON_NUM_THREADS").unwrap_or_else(|_| "mặc định".into());
    println!("[R2] Nova (nova-snark 0.76), Pallas/Vesta + IPA, RAYON_NUM_THREADS = {luong}");
    println!("[R2] hàm bước: z = (N, k) -> (N, k+1), witness (a, b) với a·b = N = 391, a, b < 2^5");

    // 1. Số ràng buộc: mạch rỗng (phần xác minh của Nova) và mạch bước của ta
    let rong = KiemPhanTich { a: 17, b: 23, rong: true, them: 0 };
    let pp0 = PublicParams::<E1, E2, KiemPhanTich>::setup(&rong, &*S1::ck_floor(), &*S2::ck_floor())
        .unwrap();
    let pp1 = PublicParams::<E1, E2, KiemPhanTich>::setup(
        &KiemPhanTich::buoc(0, 0), &*S1::ck_floor(), &*S2::ck_floor()).unwrap();
    let (c0, c1) = (pp0.num_constraints(), pp1.num_constraints());
    println!("[R2] ràng buộc (mạch chính, mạch phụ) với F rỗng: {:?}", c0);
    println!("[R2] ràng buộc (mạch chính, mạch phụ) với F = kiểm phân tích: {:?}", c1);
    println!("[R2] => F đóng góp {} ràng buộc (13 kiểm phân tích + 1 bộ đếm)", c1.0 - c0.0);

    // 2. Kiểm âm: bước 2 dùng nghiệm giả (28, 30) -- trên trường Pallas 28·30 = 840 ≠ 391
    {
        let pp = &pp1;
        let z0 = vec![F1::from(N_VD), F1::ZERO];
        let mut rs = RecursiveSNARK::new(pp, &KiemPhanTich::buoc(0, 0), &z0).unwrap();
        let mut loi_buoc = None;
        for i in 0..4 {
            let c = if i == 2 { KiemPhanTich { a: 28, b: 30, rong: false, them: 0 } } else { KiemPhanTich::buoc(i, 0) };
            if rs.prove_step(pp, &c).is_err() {
                loi_buoc = Some(i);
                break;
            }
        }
        let ok = loi_buoc.is_none() && rs.verify(pp, 4, &z0).is_ok();
        println!("[R2] kiểm âm (28, 30) ở bước 2: prove_step lỗi tại {:?}, verify = {} (mong đợi: bị từ chối)",
                 loi_buoc, ok);
        assert!(!ok);
    }

    let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/../../ket_qua/ch09");
    std::fs::create_dir_all(dir).unwrap();
    let mut csv = String::from(
        "them,so_buoc,rb_chinh,rb_phu,bien_chinh,bien_phu,setup_ms,buoc_tb_ms,buoc_tv_ms,buoc_max_ms,\
tong_ivc_ms,verify_rs_ms,nen_setup_ms,nen_prove_ms,nen_verify_ms,nen_byte,rs_byte,tai_truoc,tai_sau\n");

    // 3. Theo số bước (F = kiểm phân tích), có nén
    let ds_buoc: Vec<usize> = if nhanh { vec![10] } else { vec![10, 100, 1000] };
    let mut cau_hinh: Vec<(usize, usize, bool)> = ds_buoc.iter().map(|&n| (n, 0, true)).collect();
    // 4. Theo cỡ mạch bước (10 bước, không nén... có nén để xem cỡ proof)
    if !nhanh {
        for them in [1usize << 10, 1 << 12, 1 << 14, 1 << 16] {
            cau_hinh.push((10, them, true));
        }
    }
    for (n, them, nen) in cau_hinh {
        let k = chay(n, them, nen);
        println!(
            "[R2] them={:>6} bước={:>5} | rb={:?} | setup {:.1} ms | bước TB {:.2} ms, TV {:.2} ms, max {:.2} ms | \
IVC {:.1} ms | verify {:.2} ms | nén setup {:.1} ms, prove {:.1} ms, verify {:.2} ms, {} byte | RecursiveSNARK {} byte | tải {} -> {}",
            k.them, k.so_buoc, k.cons, k.setup_ms, k.buoc_tb, k.buoc_tv, k.buoc_max, k.tong_ivc_ms,
            k.verify_ms, k.nen_setup_ms, k.nen_prove_ms, k.nen_verify_ms, k.nen_byte, k.rs_byte,
            k.tai_truoc, k.tai_sau
        );
        csv.push_str(&format!(
            "{},{},{},{},{},{},{:.2},{:.3},{:.3},{:.3},{:.1},{:.3},{:.2},{:.2},{:.3},{},{},{},{}\n",
            k.them, k.so_buoc, k.cons.0, k.cons.1, k.vars.0, k.vars.1, k.setup_ms, k.buoc_tb,
            k.buoc_tv, k.buoc_max, k.tong_ivc_ms, k.verify_ms, k.nen_setup_ms, k.nen_prove_ms,
            k.nen_verify_ms, k.nen_byte, k.rs_byte, k.tai_truoc, k.tai_sau
        ));
        std::io::stdout().flush().ok();
    }
    let ten = if nhanh { "lab_r2_nhanh.csv" } else { "lab_r2.csv" };
    std::fs::write(format!("{dir}/{ten}"), csv).unwrap();
    println!("[R2] đã ghi ket_qua/ch09/{ten}");
}
