//! Lab R3 (bản chạy thử): MSM trên G1 và NTT trên Fr của BLS12-381, cỡ 2^10.
//! Bài lab đầy đủ (quét cỡ, so sánh luồng, đối chiếu RTL Ch.11) viết ở đợt sau.

use ark_bls12_381::{Fr, G1Affine, G1Projective};
use ark_ec::{CurveGroup, VariableBaseMSM};
use ark_ff::Zero;
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use ark_std::{
    rand::{rngs::StdRng, SeedableRng},
    UniformRand,
};
use std::time::{Duration, Instant};

const LOG_N: usize = 10;
const SO_LAN: usize = 20; // số lần đo, lấy trung vị

fn trung_vi(mut v: Vec<Duration>) -> Duration {
    v.sort();
    v[v.len() / 2]
}

fn main() {
    let n = 1usize << LOG_N;
    let mut rng = StdRng::seed_from_u64(449);
    let luong = std::thread::available_parallelism().map(|x| x.get()).unwrap_or(1);
    println!("[R3] n = 2^{} = {}, số luồng khả dụng = {}, đo {} lần, lấy trung vị", LOG_N, n, luong, SO_LAN);

    // ---- MSM: sum s_i * P_i trên G1 ----
    let diem: Vec<G1Affine> = {
        let p: Vec<G1Projective> = (0..n).map(|_| G1Projective::rand(&mut rng)).collect();
        G1Projective::normalize_batch(&p)
    };
    let vo_huong: Vec<Fr> = (0..n).map(|_| Fr::rand(&mut rng)).collect();

    let mut kq = G1Projective::zero();
    let mut t_msm = Vec::new();
    for _ in 0..SO_LAN {
        let t = Instant::now();
        kq = G1Projective::msm(&diem, &vo_huong).expect("msm lỗi");
        t_msm.push(t.elapsed());
    }
    // Đối chiếu với cách cộng ngây thơ n phép nhân vô hướng
    let t = Instant::now();
    let ngay_tho: G1Projective = diem.iter().zip(&vo_huong).map(|(p, s)| *p * s).sum();
    let t_ngay_tho = t.elapsed();
    println!("[R3] MSM G1 2^{}: {:.2?} (trung vị) | ngây thơ: {:.2?} | khớp: {}",
        LOG_N, trung_vi(t_msm), t_ngay_tho, kq == ngay_tho);

    // ---- NTT: Radix2EvaluationDomain cỡ 2^10 trên Fr ----
    let mien = Radix2EvaluationDomain::<Fr>::new(n).expect("không có miền cỡ n");
    let he_so: Vec<Fr> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
    let mut gia_tri = Vec::new();
    let mut t_fft = Vec::new();
    for _ in 0..SO_LAN {
        let t = Instant::now();
        gia_tri = mien.fft(&he_so);
        t_fft.push(t.elapsed());
    }
    let mut lai = Vec::new();
    let mut t_ifft = Vec::new();
    for _ in 0..SO_LAN {
        let t = Instant::now();
        lai = mien.ifft(&gia_tri);
        t_ifft.push(t.elapsed());
    }
    // Kiểm tra: giá trị tại điểm thứ 1 bằng đa thức tính tại omega
    let omega = mien.group_gen;
    let horner = he_so.iter().rev().fold(Fr::zero(), |acc, c| acc * omega + c);
    println!("[R3] NTT (fft) Fr 2^{}: {:.2?} | INTT (ifft): {:.2?} (trung vị)",
        LOG_N, trung_vi(t_fft), trung_vi(t_ifft));
    println!("[R3] ifft(fft(a)) == a: {} | fft(a)[1] == a(omega): {}", lai == he_so, gia_tri[1] == horner);
}
