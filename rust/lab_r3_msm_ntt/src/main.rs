//! Lab R3 (Ch.10): đo MSM trên G1 và NTT trên Fr của BLS12-381 theo cỡ n,
//! cùng thời gian một pairing (BLS12-381 và BN254, tách Miller loop / lũy thừa cuối).
//!
//! Cách chạy (script ch10/chay_r3.sh gọi hai lần: đa luồng mặc định và
//! RAYON_NUM_THREADS=1):
//!     cargo run --release -- [log_min] [log_max]
//! Mỗi dòng kết quả có dạng CSV:  CSV,phep,n,luong,trung_vi_ms,so_lan
//! Số luồng thật lấy từ rayon::current_num_threads() (tôn trọng RAYON_NUM_THREADS).

use ark_bls12_381::{Bls12_381, Fr, G1Affine, G1Projective, G2Projective};
use ark_bn254::Bn254;
use ark_ec::{pairing::Pairing, CurveGroup, PrimeGroup, VariableBaseMSM};
use ark_ff::{FftField, PrimeField, Zero};
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use ark_std::{
    rand::{rngs::StdRng, SeedableRng},
    UniformRand,
};
use std::time::Instant;

/// Trung vị (ms) của các lần đo.
fn trung_vi(mut v: Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let m = v.len();
    if m % 2 == 1 { v[m / 2] } else { (v[m / 2 - 1] + v[m / 2]) / 2.0 }
}

/// Số lần đo theo cỡ: nhỏ đo nhiều, lớn đo ít (một luồng ở cỡ >= 2^19 đo 3 lần).
fn so_lan(log_n: usize, luong: usize) -> usize {
    if log_n <= 14 { 9 } else if luong == 1 && log_n >= 19 { 3 } else { 5 }
}

fn do_ms<F: FnMut()>(lan: usize, mut f: F) -> f64 {
    let mut t = Vec::with_capacity(lan);
    for _ in 0..lan {
        let s = Instant::now();
        f();
        t.push(s.elapsed().as_secs_f64() * 1e3);
    }
    trung_vi(t)
}

fn in_tham_so() {
    // Đối chiếu với tham_so_bn_bls.py: modulus, số bit, 2-adicity lấy từ arkworks
    type Fq381 = ark_bls12_381::Fq;
    type Fq254 = ark_bn254::Fq;
    type Fr254 = ark_bn254::Fr;
    println!("[R3] BLS12-381 p = {}", Fq381::MODULUS);
    println!("[R3] BLS12-381 r = {} ({} bit, 2-adicity {})", Fr::MODULUS, Fr::MODULUS_BIT_SIZE, Fr::TWO_ADICITY);
    println!("[R3] BN254     p = {}", Fq254::MODULUS);
    println!("[R3] BN254     r = {} ({} bit, 2-adicity {})", Fr254::MODULUS, Fr254::MODULUS_BIT_SIZE, Fr254::TWO_ADICITY);
}

fn main() {
    let args: Vec<usize> = std::env::args().skip(1).filter_map(|s| s.parse().ok()).collect();
    let log_min = *args.first().unwrap_or(&10);
    let log_max = *args.get(1).unwrap_or(&20);
    let luong = rayon::current_num_threads();
    println!("[R3] n = 2^{}..2^{}, số luồng rayon = {}", log_min, log_max, luong);
    in_tham_so();

    let mut rng = StdRng::seed_from_u64(449);
    let n_max = 1usize << log_max;

    // Điểm cơ sở: cấp số cộng P_i = P0 + i*D (rẻ hơn n phép nhân vô hướng ngẫu nhiên;
    // thời gian MSM không phụ thuộc quan hệ giữa các điểm).
    let t0 = Instant::now();
    let p0 = G1Projective::rand(&mut rng);
    let d = G1Projective::rand(&mut rng);
    let mut proj = Vec::with_capacity(n_max);
    let mut cur = p0;
    for _ in 0..n_max {
        proj.push(cur);
        cur += d;
    }
    let diem: Vec<G1Affine> = G1Projective::normalize_batch(&proj);
    drop(proj);
    let vo_huong: Vec<Fr> = (0..n_max).map(|_| Fr::rand(&mut rng)).collect();
    let he_so: Vec<Fr> = (0..n_max).map(|_| Fr::rand(&mut rng)).collect();
    println!("[R3] sinh dữ liệu 2^{}: {:.1} ms", log_max, t0.elapsed().as_secs_f64() * 1e3);

    // Kiểm đúng ở cỡ nhỏ: MSM == cộng ngây thơ; ifft(fft(a)) == a; fft(a)[1] == a(omega)
    {
        let m = 1usize << log_min;
        let kq = G1Projective::msm(&diem[..m], &vo_huong[..m]).unwrap();
        let ngay_tho: G1Projective = diem[..m].iter().zip(&vo_huong[..m]).map(|(p, s)| *p * s).sum();
        let mien = Radix2EvaluationDomain::<Fr>::new(m).unwrap();
        let a = he_so[..m].to_vec();
        let v = mien.fft(&a);
        let horner = a.iter().rev().fold(Fr::zero(), |acc, c| acc * mien.group_gen + c);
        println!("[R3] kiểm 2^{}: MSM khớp ngây thơ = {}, ifft(fft(a)) == a = {}, fft(a)[1] == a(omega) = {}",
            log_min, kq == ngay_tho, mien.ifft(&v) == a, v[1] == horner);
        assert!(kq == ngay_tho && mien.ifft(&v) == a && v[1] == horner);
    }

    for log_n in log_min..=log_max {
        let n = 1usize << log_n;
        let lan = so_lan(log_n, luong);
        let ms = do_ms(lan, || {
            let r = G1Projective::msm(&diem[..n], &vo_huong[..n]).unwrap();
            let _ = std::hint::black_box(r);
        });
        println!("CSV,msm_g1,{},{},{:.4},{}", n, luong, ms, lan);

        let mien = Radix2EvaluationDomain::<Fr>::new(n).unwrap();
        let a = he_so[..n].to_vec();
        let ms = do_ms(lan, || {
            let r = mien.fft(&a);
            let _ = std::hint::black_box(r);
        });
        println!("CSV,ntt_fr,{},{},{:.4},{}", n, luong, ms, lan);
    }

    // Pairing: một pairing e(P, Q), tách Miller loop và lũy thừa cuối
    const LAN_P: usize = 51;
    let p = G1Projective::generator() * Fr::rand(&mut rng);
    let q = G2Projective::generator() * Fr::rand(&mut rng);
    let (pa, qa) = (p.into_affine(), q.into_affine());
    let ms = do_ms(LAN_P, || { let _ = std::hint::black_box(Bls12_381::pairing(pa, qa)); });
    println!("CSV,pairing_bls12_381,1,{},{:.4},{}", luong, ms, LAN_P);
    let ml = Bls12_381::multi_miller_loop([pa], [qa]);
    let ms = do_ms(LAN_P, || { let _ = std::hint::black_box(Bls12_381::multi_miller_loop([pa], [qa])); });
    println!("CSV,miller_bls12_381,1,{},{:.4},{}", luong, ms, LAN_P);
    let ms = do_ms(LAN_P, || { let _ = std::hint::black_box(Bls12_381::final_exponentiation(ml)); });
    println!("CSV,final_exp_bls12_381,1,{},{:.4},{}", luong, ms, LAN_P);
    // Song tuyến tính: e(2P, Q) == e(P, Q)^2
    let e1 = Bls12_381::pairing(pa, qa);
    let e2 = Bls12_381::pairing((p + p).into_affine(), qa);
    assert!(e2 == e1 + e1); // PairingOutput viết nhóm G_T theo lối cộng

    let p = ark_bn254::G1Projective::generator() * ark_bn254::Fr::rand(&mut rng);
    let q = ark_bn254::G2Projective::generator() * ark_bn254::Fr::rand(&mut rng);
    let (pb, qb) = (p.into_affine(), q.into_affine());
    let ms = do_ms(LAN_P, || { let _ = std::hint::black_box(Bn254::pairing(pb, qb)); });
    println!("CSV,pairing_bn254,1,{},{:.4},{}", luong, ms, LAN_P);
    let mlb = Bn254::multi_miller_loop([pb], [qb]);
    let ms = do_ms(LAN_P, || { let _ = std::hint::black_box(Bn254::multi_miller_loop([pb], [qb])); });
    println!("CSV,miller_bn254,1,{},{:.4},{}", luong, ms, LAN_P);
    let ms = do_ms(LAN_P, || { let _ = std::hint::black_box(Bn254::final_exponentiation(mlb)); });
    println!("CSV,final_exp_bn254,1,{},{:.4},{}", luong, ms, LAN_P);
    println!("[R3] song tuyến tính e(2P,Q) = e(P,Q)^2: đạt");
}
