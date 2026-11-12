"""Ch.11 (mã phụ, không in sách): hồ sơ thời gian (profile) của Prover PLONK
trong snarklib/plonk.py trên mạch N = 391 (n = 32), và ước lượng tỉ trọng MSM/NTT
cho cùng Prover ở kích thước thật từ số đo Lab R3 (ket_qua/ch10/r3.csv).

Phần 1: cProfile một lần chứng minh; tách thời gian vào ba nhóm
        MSM (opcount.msm), NTT (fft.ntt_with_root) và phần còn lại;
        đếm số MSM, tổng cỡ MSM, số NTT, tổng cỡ NTT bằng opcount.
Phần 2: ước lượng T(n) = M·t_MSM(n) + K·t_NTT(n) với M, K đếm ở phần 1 và
        t_MSM, t_NTT là trung vị đo bằng arkworks (BLS12-381, một luồng và 96 luồng).
        Đây là ƯỚC LƯỢNG: bỏ qua phép toán đa thức, băm transcript, sinh witness.
"""
import os, sys, random, cProfile, pstats, time, csv
GOC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, GOC)
from snarklib import plonk, opcount as oc
from snarklib.dsl import mach_phan_tich

c = mach_phan_tich()[0]
circ = c.compile_plonkish()
a, b, cc, x = c.plonkish_assignment({"N": 391, "a": 17, "b": 23})
pk, vk = plonk.setup(circ, tau=123)

# ---- Phần 1: profile ----
pr = cProfile.Profile()
with oc.scope() as dem:
    t0 = time.perf_counter()
    pr.enable()
    proof, _ = plonk.prove(pk, a, b, cc, x, random.Random(2019))
    pr.disable()
    tong = time.perf_counter() - t0
assert plonk.verify(vk, x, proof)
st = pstats.Stats(pr)
def cum(tenham, tenfile):
    s = 0.0
    for (f, _, fn), v in st.stats.items():
        if fn == tenham and f.endswith(tenfile):
            s += v[3]                     # thời gian tích lũy (cumtime)
    return s
t_msm = cum("msm", "opcount.py")
t_ntt = cum("ntt_with_root", "fft.py")
t_tong = sum(v[2] for v in st.stats.values())   # tổng tottime = thời gian có profile
con_lai = t_tong - t_msm - t_ntt
print("Prover PLONK, N = 391, n = %d, tau = 123, seed 2019 (Python, đường cong toy)" % pk.n)
print("số MSM = %d, tổng cỡ = %d, cỡ từng MSM = %s"
      % (dem["msm"], dem["msm_size"], dem["msm_sizes"]))
print("số phép cộng điểm (kể cả nhân đôi) = %d" % dem["add"])
print("số NTT = %d, tổng cỡ = %d" % (dem["ntt"], dem["ntt_size"]))
print("thời gian (cProfile): tổng %.3f s | MSM %.3f s (%.1f%%) | NTT %.3f s (%.1f%%) | còn lại %.3f s (%.1f%%)"
      % (t_tong, t_msm, 100 * t_msm / t_tong, t_ntt, 100 * t_ntt / t_tong,
         con_lai, 100 * con_lai / t_tong))
print("thời gian đồng hồ của lời gọi prove (có bật profile): %.3f s" % tong)

# ---- Phần 2: ước lượng ở kích thước thật từ Lab R3 ----
M, K = dem["msm"], dem["ntt"]
r3 = {}
with open(os.path.join(GOC, "ket_qua", "ch10", "r3.csv")) as f:
    for d in csv.DictReader(f):
        r3[(d["phep"], int(d["n"]), int(d["luong"]))] = float(d["trung_vi_ms"])
print("\nƯớc lượng T(n) = %d·t_MSM(n) + %d·t_NTT(n), số đo arkworks BLS12-381 (Lab R3)" % (M, K))
print("%8s %6s %10s %10s %10s %9s" % ("n", "luồng", "MSM (ms)", "NTT (ms)", "tổng (ms)", "MSM %"))
for luong in (1, 96):
    for logn in (16, 18, 20):
        n = 1 << logn
        tm, tn = r3[("msm_g1", n, luong)], r3[("ntt_fr", n, luong)]
        T = M * tm + K * tn
        print("%8s %6d %10.0f %10.0f %10.0f %8.1f%%"
              % ("2^%d" % logn, luong, M * tm, K * tn, T, 100 * M * tm / T))
