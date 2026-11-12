"""Chương 4 — Hyrax: cam kết ma trận 4 × 4 và mở MLE 16 hệ số (Mục 4.2b, bài tập lập trình).

Dùng cùng bảng z (witness N = 391) và cùng điểm r như ma_nguon_04_02_ligero.py,
nên giá trị mở phải trùng: v = z~(2, 3, 5, 7).
Nhóm: nhóm con bậc 449 của đường cong toy (chỉ phép cộng điểm).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.field import GF
from snarklib.transcript import Transcript
from snarklib.mle import eq_table, evaluate
from snarklib.pedersen import dlog_bruteforce
from snarklib.curve import toy_generator
from snarklib import hyrax, ipa

F = GF(449)
bits = lambda v: [(v >> i) & 1 for i in range(5)]
z = [1, 391, 17, 23] + bits(17) + bits(23) + [0, 0]
r = [F(2), F(3), F(5), F(7)]

pp = hyrax.setup(16)                       # 4 điểm cơ sở G_0..G_3 (độ dài một hàng)
comms, st = hyrax.commit(pp, z)            # C_i = <M[i], G>, i = 0..3
print("G =", pp.G)
print("Cam kết hàng C_0..C_3 =", comms)

L, R = eq_table(F, r[:2]), eq_table(F, r[2:])
w = [sum((L[i] * st["M"][i][j] for i in range(4)), F(0)) for j in range(4)]
print("L = eq(·, (2,3)) =", L, "  R = eq(·, (5,7)) =", R)
print("w = L^T·M =", w, "  v = <w, R> =", ipa.inner(w, R))
C_w = hyrax.combine_rows(comms, L)         # Verifier tự tính từ 4 cam kết
print("C_w = Σ L_i·C_i =", C_w)

v, proof = hyrax.open(pp, st, r, Transcript(b"ch4-hyrax"))
cnt = ipa.VCount()
ok = hyrax.verify(pp, comms, r, v, proof, Transcript(b"ch4-hyrax"), 16, cnt=cnt)
print("v =", v, " (mle.evaluate:", evaluate(F, z, r), ")  Verifier chấp nhận:", ok)
print("Kích thước:", hyrax.sizes(16), " Verifier:", cnt.as_dict())
bad = hyrax.verify(pp, comms, r, v + 1, proof, Transcript(b"ch4-hyrax"), 16)
print("Khai sai v + 1 -> chấp nhận:", bad)

# Cảnh báo: nhóm toy chỉ có 449 phần tử, log rời rạc của G_0 theo G1 tìm ngay
_, _, G1 = toy_generator()
print("log_{G1}(G_0) =", dlog_bruteforce(G1, pp.G[0]), " (vét cạn — toy không có binding thật)")
assert ok and not bad and v == evaluate(F, z, r)
