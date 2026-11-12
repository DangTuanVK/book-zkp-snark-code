"""Kiểm các tham số xuyên suốt (Chương 1): r, g, ω_16, #E, G1, bậc nhúng và
tính song tuyến tính e(aG, bG) = e(G, G)^(ab) của pairing toy."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib import GF, is_prime, factorize
from snarklib.curve import toy_generator
from snarklib.pairing import e, embedding_degree

# 1. Trường vô hướng F_r
r = 449
F = GF(r)
assert is_prime(r)
print("r = %d nguyên tố: %s, r - 1 = %s" % (r, is_prime(r), factorize(r - 1)))
g = F.generator()
w16 = F.root_of_unity(16)                 # = g^((r-1)/16) = g^28
assert int(g) == 3 and int(w16) == 349
assert w16 ** 16 == 1 and w16 ** 8 != 1   # bậc đúng bằng 16
print("g = %d, ω_16 = g^28 = %d, bậc ω_16 = %d"
      % (int(g), int(w16), w16.order()))

# 2. Đường cong toy E: y^2 = x^3 + x trên F_q
E, P, G = toy_generator()
q = E.q
nE = E.count_points()                     # duyệt mọi x trong F_q
assert nE == q + 1 == 12 * r              # siêu kỳ dị: #E = q + 1
print("%s, #E = %d = %s" % (E, nE, factorize(nE)))
print("P = (%d, %d), G1 = 12·P = (%d, %d)"
      % (int(P.x), int(P.y), int(G.x), int(G.y)))
assert not G.is_infinity() and (r * G).is_infinity()
k = embedding_degree(q, r)
assert k == 2
print("bậc G1 = %d, bậc nhúng k = %d" % (G.order(nE), k))

# 3. Pairing toy: không suy biến và song tuyến tính
gT = e(G, G)
assert gT != 1 and gT ** r == 1           # e(G, G) có bậc r trong F_{q^2}^*
print("e(G1, G1) = %s, e(G1, G1)^r = %s" % (gT, gT ** r))
for a, b in [(2, 3), (17, 23), (28, 30), (448, 448)]:
    trai = e(a * G, b * G)
    phai = gT ** (a * b % r)
    assert trai == phai
    print("a = %3d, b = %3d: e(aG, bG) = %-16s = e(G,G)^%d"
          % (a, b, trai, a * b % r))
print("Mọi kiểm tra đạt.")
