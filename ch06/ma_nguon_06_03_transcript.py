"""Chương 6 -- Transcript thật của mini-PLONK cho N = 391 (tab:ch6-cac-vong).

In từng vòng: cam kết, thách thức β, γ, α, ζ, v, u, giá trị mở; z(ω^i);
kích thước proof; số phép toán Prover/Verifier đo bằng snarklib.opcount.
SRS: τ = 123 (chỉ để minh họa), bậc n + 5 = 37. Seed làm mù cố định.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib import plonk, opcount as oc
from snarklib.dsl import mach_phan_tich

c = mach_phan_tich()[0]
circ = c.compile_plonkish()
a, b, cc, x = c.plonkish_assignment({"N": 391, "a": 17, "b": 23})
pk, vk = plonk.setup(circ, tau=123)
n = pk.n
print("n = %d, ω = %d, SRS: %d điểm G1 ([τ^0..τ^%d]_1) + 2 điểm G2"
      % (n, int(pk.omega), len(pk.srs.g1), len(pk.srs.g1) - 1))

with oc.scope() as cp:
    P = plonk.Prover(pk, a, b, cc, x, random.Random(2019))
    proof = P.run()
print("Prover dừng-làm-lại vì mẫu số z bằng 0:", P.aborts, "lần")
vong = 0
for rnd, ten, val in P.log:
    if rnd != vong:
        vong = rnd
        print("--- Vòng %d" % rnd)
    print("  %-16s %s" % (ten, val))

print("z(ω^i), i = 0..%d:" % (n - 1))
print(" ", [int(v) for v in P.z_vals])
print("z(ω^20..31) =", sorted(set(int(v) for v in P.z_vals[20:])),
      "(mọi vòng σ nằm trong hàng 0..19 nên tích đã quay về 1 tại hàng 20)")
print("bậc a, b, c, z, t_lo, t_mid, t_hi:", P.a_p.degree(), P.b_p.degree(), P.c_p.degree(),
      P.z_p.degree(), P.t_lo.degree(), P.t_mid.degree(), P.t_hi.degree())

with oc.scope() as cv:
    ok = plonk.verify(vk, x, proof)
print("verify(N = 391):", ok)
print("kích thước proof: %d phần tử G1 + %d phần tử F_r" % (len(plonk.G1_FIELDS), len(plonk.FR_FIELDS)))
for ten_ct, g1, fr in (("BN254", 32, 32), ("BLS12-381", 48, 32)):
    print("  nếu nén điểm trên %s: 9·%d + 6·%d = %d byte" % (ten_ct, g1, fr, 9 * g1 + 6 * fr))
print("Prover (một lần chạy thành công + %d lần làm lại):" % P.aborts)
print("  MSM: %d lần, cỡ %s, tổng cỡ %d" % (cp["msm"], cp["msm_sizes"], cp["msm_size"]))
print("  NTT cỡ %d: %d lần; pairing: %d" % (n, cp.get("ntt", 0), cp.get("pairing", 0)))
print("Verifier: MSM cỡ %s, vô hướng khác 0/1 = %d, pairing = %d, lũy thừa cuối = %d"
      % (cv["msm_sizes"], cv["msm_nontrivial"], cv["pairing"], cv["final_exp"]))
assert ok and P.z_vals[0] == 1 and {int(v) for v in P.z_vals[20:]} == {1} and cv["pairing"] == 2 and cv["msm_nontrivial"] == 18
assert not plonk.verify(vk, [389], proof)
