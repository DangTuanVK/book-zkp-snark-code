"""Chương 9 — Mã phụ cho Mục 9.2c (không in trong sách): extractor của folding Nova.

(1) Extractor theo Bổ đề 7 của bài Nova: từ 3 transcript chung số hạng chéo T,
    khác thách thức r, nội suy (W1, W2) từ 2 transcript và (E1, T, E2) từ 3
    transcript; kiểm hai witness trích được thỏa hai instance ban đầu.
(2) Sai số soundness 2/|F|: instance thứ hai sai (17, 24); với mỗi giá trị v
    Prover khai thay cho T ở hàng 12, đếm số thách thức r in F_449 làm instance
    gộp (với witness gộp trung thực) thỏa R1CS nới lỏng.
(3) Mài Fiat-Shamir trên trường toy: thử v_12 cho tới khi thách thức băm ra
    đúng giá trị r* mà v_12 cần.
Không cam kết (witness rõ), z = (u, x, W) theo snarklib.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.dsl import mach_phan_tich
from snarklib.folding import RelaxedR1CS
from snarklib.transcript import Transcript
from snarklib import sumcheck as sc

c = mach_phan_tich()[0]
R = RelaxedR1CS(c.compile_r1cs())
F = R.F
wit = lambda a, b: R.trivial(*c.r1cs_assignment({"N": 391, "a": a, "b": b}))

# Instance thứ nhất: instance nới lỏng của Ví dụ ex:ch9-nova-391 (u = 398, E != 0)
(i1, W1), (i2, W2) = wit(17, 23), wit(23, 17)
U, WU, _, r0 = R.fold(i1, W1, i2, W2, Transcript(b"ch9-nova"))
u_new, W_new = wit(17, 23)                        # instance mới, u = 1, E = 0
print("(1) Gộp U (u = %d, E != 0) với instance mới (17, 23)" % U["u"])
T = R.cross_term(U, WU, u_new, W_new)

def gap(r):                                        # Prover trung thực với thách thức r
    r = F(r)
    inst = {"u": U["u"] + r * u_new["u"], "x": [a + r * b for a, b in zip(U["x"], u_new["x"])],
            "E": [e1 + r * t + r * r * e2 for e1, t, e2 in zip(U["E"], T, u_new["E"])]}
    W = [a + r * b for a, b in zip(WU, W_new)]
    assert R.is_satisfied(inst, W)
    return inst, W

rs = [F(3), F(10), F(100)]
tr = [gap(r) for r in rs]
# W: nội suy bậc 1 từ hai transcript đầu
ra, rb = rs[0], rs[1]
W2x = [(wa - wb) / (ra - rb) for wa, wb in zip(tr[0][1], tr[1][1])]
W1x = [wa - ra * w2 for wa, w2 in zip(tr[0][1], W2x)]
# E: nội suy bậc 2 (Vandermonde 3x3) từ ba transcript
def noi_suy_bac2(ys):
    # tìm (c0, c1, c2) với c0 + c1 r_i + c2 r_i^2 = y_i, Lagrange
    out = [F(0)] * 3
    for i in range(3):
        others = [rs[j] for j in range(3) if j != i]
        den = (rs[i] - others[0]) * (rs[i] - others[1])
        # (X - o0)(X - o1) = X^2 - (o0 + o1) X + o0 o1
        co = [others[0] * others[1], -(others[0] + others[1]), F(1)]
        for k in range(3):
            out[k] = out[k] + ys[i] * co[k] / den
    return out
E1x, Tx, E2x = [], [], []
for j in range(len(T)):
    c0, c1, c2 = noi_suy_bac2([tr[i][0]["E"][j] for i in range(3)])
    E1x.append(c0); Tx.append(c1); E2x.append(c2)
print("    thách thức r =", [int(r) for r in rs])
print("    W2 trích (a, b) =", [int(v) for v in W2x[:2]], "| W1 trích (a, b) =", [int(v) for v in W1x[:2]])
print("    E1 trích =", [int(v) for v in E1x])
print("    E2 trích khác 0:", sum(1 for v in E2x if v), "| T trích khớp T:", Tx == T)
ok1 = R.is_satisfied({"u": U["u"], "x": U["x"], "E": E1x}, W1x) and W1x == WU and E1x == U["E"]
ok2 = R.is_satisfied({"u": u_new["u"], "x": u_new["x"], "E": E2x}, W2x) and W2x == W_new
print("    witness trích được thỏa U:", ok1, "| thỏa instance mới:", ok2)
assert ok1 and ok2
# Hai transcript không đủ cho E: mọi (E1, T, E2) dạng (c0 + k r_a r_b, c1 - k(r_a + r_b), c2 + k) đều khớp
print("    hai transcript: họ nghiệm E có 449 phần tử (một tham số tự do k in F_449)")

# (2) Instance thứ hai sai: (17, 24), 17·24 - 391 = 17
ib, Wb = wit(17, 24)
Tb = R.cross_term(i1, W1, ib, Wb)
res = R.residual(ib, Wb)
print("(2) Instance sai (17, 24): phần dư khác 0 ở hàng", [j for j, v in enumerate(res) if v],
      "giá trị", [int(v) for v in res if v])
def so_r_dat(v12):
    v = list(Tb); v[12] = F(v12)
    dem = []
    for r in range(449):
        rr = F(r)
        inst = {"u": i1["u"] + rr * ib["u"], "x": [a + rr * b for a, b in zip(i1["x"], ib["x"])],
                "E": [e1 + rr * t + rr * rr * e2 for e1, t, e2 in zip(i1["E"], v, ib["E"])]}
        W = [a + rr * b for a, b in zip(W1, Wb)]
        if R.is_satisfied(inst, W):
            dem.append(r)
    return dem
hist = {}
for v12 in range(449):
    k = len(so_r_dat(v12))
    hist[k] = hist.get(k, 0) + 1
print("    số v_12 theo số thách thức r đạt:", dict(sorted(hist.items())))
print("    v_12 = T_12 =", int(Tb[12]), "-> r đạt:", so_r_dat(int(Tb[12])),
      "| v_12 = 0 -> r đạt:", so_r_dat(0))
assert max(hist) <= 2

# (3) Mài Fiat-Shamir: r* = (v_12 - T_12)/17 phải trùng thách thức băm
lan, trung = 0, []
for v12 in range(449):
    lan += 1
    v = list(Tb); v[12] = F(v12)
    t = Transcript(b"ch9-mai")
    for nm, inst in ((b"/U1", i1), (b"/U2", ib)):
        t.absorb(b"nova" + nm, [int(inst["u"])] + [int(x) for x in inst["x"]] + [int(x) for x in inst["E"]])
    t.absorb(b"nova/T", [int(x) for x in v])
    r = sc._challenge(F, t, b"nova/r")
    if r != 0 and (F(v12) - Tb[12]) == F(17) * r:
        trung.append((v12, int(r)))
print("(3) Mài Fiat-Shamir qua %d giá trị v_12: trúng %d lần %s" % (lan, len(trung), trung))
