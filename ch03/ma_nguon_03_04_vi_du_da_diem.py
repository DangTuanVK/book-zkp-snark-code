"""Ví dụ số Mục 3.2 (một đa thức, nhiều điểm): mở một đa thức bậc 3 tại hai điểm trên F_449.

f(X) = X^3 + 17X + 23, S = {1, 2}, cửa sập τ = 10 (chỉ để minh họa).
In Z_S, r(X), thương ψ(X) = (f - r)/Z_S và kiểm phương trình pairing
    e(C - [r(τ)]_1, [1]_2) = e(W, [Z_S(τ)]_2).
Đường cong toy là Type 1 (G2 = G1).
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from snarklib import kzg
from snarklib.poly import Poly, lagrange_interpolate
from snarklib.pairing import e

F = kzg.F
TAU = 10
f = Poly(F, [23, 17, 0, 1])                  # 23 + 17X + X^3
S = [1, 2]

# --- phía prover: Z_S, r(X), thương ψ(X) -------------------------------------
ys = [f(F(s)) for s in S]                     # f(1), f(2)
ZS = kzg.vanishing(S)                         # (X - 1)(X - 2)
r = lagrange_interpolate(F, S, ys)            # đa thức bậc < 2 đi qua (s, f(s))
psi, du = (f - r).divmod(ZS)
assert du.is_zero()

# --- cùng ví dụ, mở một điểm z = 1 để so sánh ---------------------------------
q1, du1 = (f - ys[0]).divmod(Poly(F, [-1, 1]))
assert du1.is_zero()

# --- cam kết trên đường cong toy với τ = 10 ----------------------------------
srs = kzg.setup(3, TAU, g2_degree=2)
C = kzg.commit(srs, f)
ys2, W = kzg.open_multi(srs, f, S)
assert ys2 == ys
G = srs.G1
t = F(TAU)

print("f(X)        =", f)
print("S           =", S, " f(S) =", [int(y) for y in ys])
print("Z_S(X)      =", ZS)
print("r(X)        =", r)
print("psi(X)      =", psi, "  (thương (f - r)/Z_S)")
print("q_1(X)      =", q1, "  (thương mở một điểm z = 1)")
print("tau = %d: f(tau) = %d, r(tau) = %d, Z_S(tau) = %d, psi(tau) = %d"
      % (TAU, int(f(t)), int(r(t)), int(ZS(t)), int(psi(t))))
print("kiểm số: psi(tau)·Z_S(tau) + r(tau) = %d = f(tau)" % int(psi(t) * ZS(t) + r(t)))
print("G1          =", G)
print("C = [f(tau)]_1 =", C)
print("W = [psi(tau)]_1 =", W)
print("[r(tau)]_1  =", kzg.commit(srs, r))
print("[Z_S(tau)]_2 =", int(ZS(t)) * srs.G2)
lhs = e(C - kzg.commit(srs, r), srs.G2)
rhs = e(W, int(ZS(t)) * srs.G2)
print("e(C - [r(tau)]_1, [1]_2) =", lhs)
print("e(W, [Z_S(tau)]_2)       =", rhs)
ok = kzg.verify_multi(srs, C, S, ys, W)
bad = kzg.verify_multi(srs, C, S, [ys[0], ys[1] + 1], W)
print("verify_multi(đúng) =", ok, "; verify_multi(f(2) sai lệch 1) =", bad)
assert lhs == rhs and ok and not bad
assert f(t) == psi(t) * ZS(t) + r(t)
