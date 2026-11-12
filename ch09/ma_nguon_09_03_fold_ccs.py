"""Chương 9 — Một bước multi-folding kiểu HyperNova trên mạch ví dụ (lst:ch9-fold-ccs).

LCCCS tích lũy (witness (17, 23)) + CCCS mới (witness (23, 17)) -> LCCCS mới, gộp
bằng một sum-check bậc d + 1 = 3 trên s = 4 biến. Đơn giản hóa so với HyperNova:
không cam kết (witness rõ), không ZK; z = (u, x, w) theo thứ tự của snarklib.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.dsl import mach_phan_tich
from snarklib.folding import CCSPadded
from snarklib.transcript import Transcript

c = mach_phan_tich()[0]
P = CCSPadded(c.compile_ccs())             # CCS bậc 2 từ R1CS: M = (A, B, C), S = ({0,1}, {2})
print("CCS: m = %d -> %d hàng, n = %d -> %d cột, s = %d, t = %d, bậc d = %d"
      % (P.ccs.m, P.m, P.ccs.n, P.n, P.s, P.t, P.d))
wit = lambda a, b: c.r1cs_assignment({"N": 391, "a": a, "b": b})
x1, w1 = wit(17, 23)
x2, w2 = wit(23, 17)

# Prover: tuyến tính hóa instance đầu, sau đó gộp với CCCS mới
tp = Transcript(b"ch9-hypernova")
L1 = P.linearize(x1, w1, tp)
print("LCCCS_1: u =", L1["u"], " x =", L1["x"], " r =", L1["r"], " v =", L1["v"])
proof, L, w_new = P.hypernova_prove(L1, w1, x2, w2, tp)
for j, s in enumerate(proof["sc"]["polys"]):
    print("  vòng %d: s_%d(0..3) =" % (j + 1, j + 1), s)
print("sigma =", proof["sigma"], " theta =", proof["theta"])
print("LCCCS mới: u =", L["u"], " x =", L["x"], " r' =", L["r"], " v' =", L["v"])
print("w' (a, b) =", w_new[:2], "= w1 + rho·w2")
assert P.lcccs_is_satisfied(L, w_new)

# Verifier: chỉ thấy instance + proof, tự suy ra cùng LCCCS mới
tv = Transcript(b"ch9-hypernova")
ok, Lv = P.hypernova_verify(P.linearize(x1, w1, tv), x2, proof, tv)
print("Verifier chấp nhận:", ok, "| LCCCS khớp Prover:", Lv == L)
assert ok and Lv == L
n_sc = sum(len(s) for s in proof["sc"]["polys"])
print("Kích thước proof gộp: %d (sum-check) + %d (sigma, theta) = %d phần tử F_449; bỏ s_j(1): %d"
      % (n_sc, 2 * P.t, n_sc + 2 * P.t, n_sc - P.s + 2 * P.t))

# Gộp tiếp bước 2 (LCCCS vừa có + CCCS (17, 23)) để thấy chuỗi IVC
proof2, L2, w2n = P.hypernova_prove(L, w_new, x1, w1, tp, label=b"hn2")
print("Bước 2: u =", L2["u"], " thỏa:", P.lcccs_is_satisfied(L2, w2n))
assert P.lcccs_is_satisfied(L2, w2n)

# Kiểm âm: CCCS sai (a·b = 17·24 ≠ 391) bị loại; nghiệm giả (28, 30) lọt (Ch.12)
def thu(a, b, k=20):
    xb, wb = wit(a, b)
    dat = 0
    for i in range(k):
        lab = b"neg%d" % i
        tp, tv = Transcript(lab), Transcript(lab)
        pf, Lb, wn = P.hypernova_prove(P.linearize(x1, w1, tp), w1, xb, wb, tp)
        okb, _ = P.hypernova_verify(P.linearize(x1, w1, tv), xb, pf, tv)
        dat += okb and P.lcccs_is_satisfied(Lb, wn)
    return dat
sai, gia = thu(17, 24), thu(28, 30)
print("Witness (17, 24): đạt %d/20 transcript | nghiệm giả (28, 30): đạt %d/20" % (sai, gia))
assert sai <= 2 and gia == 20
