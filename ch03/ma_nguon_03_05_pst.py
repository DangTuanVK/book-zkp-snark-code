"""KZG đa tuyến tính (PST13) — ví dụ 2 biến tính tay được và kiểm 3 biến (Mục 3.2b).

Bảng 2 biến (big-endian, x_1 là bit cao): q(0,0)=17, q(0,1)=23, q(1,0)=2, q(1,1)=31.
Cửa sập τ = (5, 6), điểm mở z = (3, 4). Toy: Type 1 (G2 = G1).
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from snarklib import pst, opcount as oc
from snarklib.mle import evaluate, eq_table

F = pst.F
T = [17, 23, 2, 31]
TAU, Z = [5, 6], [3, 4]

srs = pst.setup(2, TAU)
v = evaluate(F, T, Z)
w1, w2 = pst.quotients(T, Z)          # w_1(X_2) và hằng số w_2
print("bảng q             =", T, " q(z) = q(3, 4) =", int(v))
print("w_1 trên X_2∈{0,1} =", [int(x) for x in w1], " (= q(1, X_2) - q(0, X_2))")
print("w_2 (hằng)         =", [int(x) for x in w2], " (= hiệu sau khi cố định x_1 = 3)")
# kiểm Fact 15.1 tại một điểm X tùy ý: q(X) - v = (X1 - z1) w1(X2) + (X2 - z2) w2
X = [F(100), F(200)]
lhs = evaluate(F, T, X) - v
rhs = (X[0] - Z[0]) * evaluate(F, w1, X[1:]) + (X[1] - Z[1]) * w2[0]
print("kiểm q(X) - v = Σ (X_i - z_i) w_i tại X = (100, 200):", int(lhs), "=", int(rhs))
print("SRS tables[0] = [eq(b, tau)]_1 với eq(b, tau) =", [int(x) for x in eq_table(F, TAU)])

C = pst.commit(srs, T)
with oc.scope() as cp:
    v2, Ws = pst.open(srs, T, Z)
with oc.scope() as cv:
    ok = pst.verify(srs, C, Z, v2, Ws)
print("C =", C, " W_1 =", Ws[0], " W_2 =", Ws[1])
print("verify =", ok, "; prover MSM cỡ", cp["msm_sizes"], "; verifier", oc.summary(cv))
print("verify với v sai lệch 1 =", pst.verify(srs, C, Z, v2 + 1, Ws))

# 3 biến: so với snarklib.mle trên mọi điểm Boolean và một điểm ngẫu nhiên
T3, Z3 = [17, 23, 391, 5, 0, 1, 2, 448], [10, 20, 30]
srs3 = pst.setup(3, [7, 8, 9])
C3 = pst.commit(srs3, T3)
v3, W3 = pst.open(srs3, T3, Z3)
print("3 biến: q(10,20,30) =", int(v3), "; proof", len(W3), "phần tử G1; verify =",
      pst.verify(srs3, C3, Z3, v3, W3))
assert ok and lhs == rhs and v3 == evaluate(F, T3, Z3) and cv["pairing"] == 3
assert not pst.verify(srs, C, Z, v2 + 1, Ws)
