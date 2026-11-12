"""Chương 6 -- Kiểm chứng phụ cho mục Marlin (mã hóa row/col/val), không in trong sách.

R1CS 13 ràng buộc x 14 biến của mạch N = 391 (Ch.5). Theo chiesa2020marlin §2.4:
H đánh chỉ số hàng/cột (|H| = 16, ω_16 = 349), K đánh chỉ số các ô khác 0.
  M̂(X, Y) = Σ_{κ∈K} u_H(X, row(κ)) · u_H(Y, col(κ)) · val(κ),
  u_H(X, Y) = (v_H(X) - v_H(Y)) / (X - Y), u_H(a, a) = n·a^{n-1},
và val(κ) = M[r][c] / (u_H(r, r) · u_H(c, c)) (chú thích 10 của bài).
Kiểm M̂(h1, h2) = M[h1][h2] trên toàn bộ H × H, đếm số đa thức indexer xuất ra.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.field import GF
from snarklib.dsl import mach_phan_tich

F = GF(449)
c = mach_phan_tich()[0]
r1 = c.compile_r1cs()
nH = 16
om = F(3) ** (448 // nH)
H = [om ** i for i in range(nH)]


def uH(X, Y):
    if X == Y:
        return F(nH) * X ** (nH - 1)
    return (X ** nH - Y ** nH) / (X - Y)


nnz = [sum(1 for row in M for v in row if v) for M in (r1.A, r1.B, r1.C)]
nK = 1
while nK < max(nnz):
    nK *= 2
print("R1CS: %d ràng buộc x %d biến; ω_16 = %d; số ô khác 0 của A, B, C = %s -> |K| = %d"
      % (r1.m, r1.n, int(om), nnz, nK))
for ten, M in (("A", r1.A), ("B", r1.B), ("C", r1.C)):
    o = [(i, j, F(v)) for i, row in enumerate(M) for j, v in enumerate(row) if v]
    row = [H[i] for i, _, _ in o]
    col = [H[j] for _, j, _ in o]
    val = [v / (uH(H[i], H[i]) * uH(H[j], H[j])) for i, j, v in o]
    dung = all(
        sum((uH(H[h1], rw) * uH(H[h2], cl) * vl for rw, cl, vl in zip(row, col, val)), F(0))
        == (F(M[h1][h2]) if h1 < r1.m and h2 < r1.n else F(0))
        for h1 in range(nH) for h2 in range(nH))
    print("  %s: %2d ô; M̂ khôi phục đúng mọi ô trên H x H: %s" % (ten, len(o), dung))
    if ten == "A":
        i, j, v = o[-1]
        print("     ví dụ ô cuối của A: hàng %d, cột %d, giá trị %d -> row = ω^%d = %d, col = ω^%d = %d,"
              " u_H(r,r) = %d, u_H(c,c) = %d, val = %d"
              % (i, j, int(v), i, int(H[i]), j, int(H[j]), int(uH(H[i], H[i])),
                 int(uH(H[j], H[j])), int(val[-1])))
    assert dung
print("indexer Marlin: 3 ma trận x 3 đa thức (row, col, val) bậc < |K| = %d -> 9 đa thức" % nK)
print("so với PLONK (cùng mạch): 8 đa thức tiền xử lý bậc < n = 32")
print("số ô khác 0 tổng cộng: %d; số ô của ma trận đầy đủ 16 x 16 x 3 = %d" % (sum(nnz), 3 * nH * nH))
