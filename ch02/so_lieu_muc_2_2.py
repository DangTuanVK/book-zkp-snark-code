"""Chương 2, Mục 2.2 -- số liệu cho các chứng minh mẫu của KZG (tệp phụ, không in).

(1) Phá tính ràng buộc theo giá trị đánh giá trên SRS toy (tau = 10 đã biết):
    hai lời mở khác nhau tại cùng một điểm z cho cùng cam kết C.
(2) Phép quy dẫn của Định lý ràng buộc: từ hai lời mở suy ra lời giải q-SDH.
(3) Đối thủ đại số trong AGM: phép kiểm lý tưởng sai, phép kiểm thật đúng,
    phân tích đa thức hiệu cho ra tau (phép quy dẫn về q-DLOG).
(4) Lời mở gộp nhiều đa thức: tỉ lệ gamma làm tổng sai lệch triệt tiêu.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.field import GF
from snarklib.poly import Poly
from snarklib.pairing import e
from snarklib import kzg

r = 449
F = GF(r)
tau = 10                                   # cửa sập toy, như Ví dụ 2.x
srs = kzg.setup(3, tau)
G = srs.G1
f = Poly(F, [23, 17, 0, 1])                # f(X) = X^3 + 17X + 23
C = kzg.commit(srs, f)
print("== (1) SRS bậc 3, tau = %d; G1 = %s; C = %s" % (tau, G, C))
assert str(C) == str(G * 295)

z = 3
y, W = kzg.open(srs, f, z)                 # lời mở trung thực
print("lời mở trung thực tại z = %d: y = %d, W = %s, verify = %s"
      % (z, int(y), W, kzg.verify(srs, C, z, y, W)))
y2 = y + 1                                 # giá trị giả
# đối thủ biết tau: W' = [(f(tau) - y') / (tau - z)]_1
w2 = (F(295) - y2) / (F(tau) - F(z))
W2 = G * int(w2)
print("lời mở giả: y' = %d, W' = %d·G1 = %s, verify = %s"
      % (int(y2), int(w2), W2, kzg.verify(srs, C, z, y2, W2)))
assert kzg.verify(srs, C, z, y, W) and kzg.verify(srs, C, z, y2, W2)
wq = (F(295) - y) / (F(tau) - F(z))
print("W = %d·G1 (thương q(tau) = %d)" % (int(wq), int(wq)))

print("\n== (2) Quy dẫn: P = (y' - y)^(-1) · (W - W'), c = -z")
d = int(y2 - y)
P = (W - W2) * int(F(d).inverse())
c = int(-F(z))
print("c = %d, P = %s" % (c, P))
inv = F(tau + c).inverse()
print("1/(tau + c) mod 449 = %d; P == [1/(tau+c)]_1: %s"
      % (int(inv), P == G * int(inv)))
lhs = e(P, srs.g2[1] + srs.G2 * c)
rhs = e(G, srs.G2)
print("e(P, [tau]_2 + c[1]_2) = %s; e(G1, G1) = %s" % (lhs, rhs))
assert lhs == rhs

print("\n== (3) AGM: biểu diễn của W' là hằng số w' = %d" % int(w2))
# phép kiểm lý tưởng: f(X) - y' == H(X)(X - z) với H(X) = w' (hằng)
Rpoly = f - Poly(F, [y2]) - Poly(F, [w2]) * Poly(F, [-F(z), 1])
print("R(X) = f(X) - y' - w'(X - z) =", Rpoly)
roots = [t for t in range(r) if Rpoly(F(t)) == 0]
print("nghiệm của R trong F_449:", roots)
cand = [t for t in roots if G * t == srs.g1[1]]
print("nghiệm t có [t]_1 = [tau]_1:", cand)
assert cand == [tau]
# lời mở trung thực: phép kiểm lý tưởng đúng
Rh = f - Poly(F, [y]) - (f - Poly(F, [y])) // Poly(F, [-F(z), 1]) * Poly(F, [-F(z), 1])
print("lời mở trung thực: hiệu lý tưởng bằng 0:", Rh.is_zero())

print("\n== (4) Gộp t = 2 đa thức tại cùng z: tổng sai lệch triệt tiêu theo gamma")
g = Poly(F, [5, 0, 1])                     # g(X) = X^2 + 5
s1, s2 = f(F(z)) + 1, g(F(z)) + 7          # hai giá trị sai
bad = [gm for gm in range(r) if (f(F(z)) - s1) + F(gm) * (g(F(z)) - s2) == 0]
print("số gamma làm (f(z)-s1) + gamma(g(z)-s2) = 0:", len(bad), bad)
assert len(bad) <= 1
