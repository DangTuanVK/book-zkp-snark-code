"""Chương 6 -- PLONK và Groth16 trên cùng mạch N = 391 (tab:ch6-plonk-vs-groth16).

PLONK: đếm bằng snarklib.opcount. Groth16 (snarklib.groth16, Ch.5) không đi qua
opcount nên ta bọc tạm hàm msm, e và phép nhân điểm của module để đếm.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib import plonk, groth16, opcount as oc
from snarklib.curve import Point, toy_generator
from snarklib.dsl import mach_phan_tich

c = mach_phan_tich()[0]
inp = {"N": 391, "a": 17, "b": 23}
# --- PLONK
circ = c.compile_plonkish()
a, b, cc, x = c.plonkish_assignment(inp)
pk, vk = plonk.setup(circ, tau=123)
with oc.scope() as pp:
    proof, _ = plonk.prove(pk, a, b, cc, x, random.Random(2019))
with oc.scope() as pv:
    assert plonk.verify(vk, x, proof)
# --- Groth16: bọc msm, e, Point.__mul__ để đếm
dem = {"msm": [], "e": 0, "smul": 0}
goc_msm, goc_e, goc_mul = groth16.msm, groth16.e, Point.__mul__
def msm_dem(s, P):
    dem["msm"].append(len(s)); return goc_msm(s, P)
def e_dem(P, Q):
    dem["e"] += 1; return goc_e(P, Q)
def mul_dem(P, k):
    dem["smul"] += 1; return goc_mul(P, k)
G = toy_generator()[2]
r1 = c.compile_r1cs()
qap = groth16.r1cs_to_qap(r1)
xg, w = c.r1cs_assignment(inp)
gpk, gvk, _ = groth16.setup(qap, G, G, random.Random(2016))
groth16.msm, groth16.e, Point.__mul__ = msm_dem, e_dem, mul_dem
Point.__rmul__ = mul_dem
try:
    gp = groth16.prove(gpk, qap, [1] + xg + w, random.Random(5))
    g_prove = dict(msm=list(dem["msm"]), smul=dem["smul"])
    dem.update(msm=[], e=0, smul=0)
    assert groth16.verify(gvk, xg, gp)
    g_ver = dict(dem)
finally:
    groth16.msm, groth16.e, Point.__mul__, Point.__rmul__ = goc_msm, goc_e, goc_mul, goc_mul

print("%-34s %-26s %s" % ("", "Groth16", "PLONK"))
print("%-34s %-26s %s" % ("số hàng / miền", "R1CS 13 (+2) -> n = %d" % qap.n,
                          "Plonkish 20 -> n = %d" % pk.n))
print("%-34s %-26s %s" % ("setup", "riêng cho mạch", "SRS phổ quát bậc %d" % (pk.n + 5)))
print("%-34s %-26s %s" % ("proof", "2 G1 + 1 G2", "9 G1 + 6 F_r"))
print("%-34s %-26s %s" % ("byte (nén, BLS12-381)", "2·48 + 96 = 192", "9·48 + 6·32 = 624"))
print("%-34s %-26s %s" % ("byte (nén, BN254)", "2·32 + 64 = 128", "9·32 + 6·32 = 480"))
print("%-34s %-26s %s" % ("Prover: số MSM (cỡ)", "%d %s" % (len(g_prove["msm"]), g_prove["msm"]),
                          "%d %s" % (pp["msm"], pp["msm_sizes"])))
print("%-34s %-26s %s" % ("Prover: tổng cỡ MSM", sum(g_prove["msm"]), pp["msm_size"]))
print("%-34s %-26s %s" % ("Prover: nhân vô hướng (mọi loại)", g_prove["smul"],
                          "%d (trong MSM)" % pp["msm_nontrivial"]))
print("%-34s %-26s %s" % ("Prover: NTT cỡ n", "(INTT trong QAP)", "%d" % pp["ntt"]))
print("%-34s %-26s %s" % ("Verifier: pairing", "%d (+ e(α,β) tính sẵn)" % g_ver["e"],
                          "%d" % pv["pairing"]))
print("%-34s %-26s %s" % ("Verifier: MSM (cỡ)", g_ver["msm"], pv["msm_sizes"]))
print("%-34s %-26s %s" % ("Verifier: nhân vô hướng", g_ver["smul"], pv["msm_nontrivial"]))
assert g_ver["e"] == 3 and pv["pairing"] == 2 and len(plonk.Proof._fields) == 15
