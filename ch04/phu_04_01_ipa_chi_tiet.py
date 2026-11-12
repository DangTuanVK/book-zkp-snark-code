"""Chương 4 — Số liệu chi tiết cho Mục 4.2a/4.2b (không in trong sách).

(1) Kiểm tay một vòng gập của ma_nguon_04_01: x, x^-1, x^2, x^-2, P, P'.
(2) PCS-IPA n = 4 (transcript 'ch4-pcs'): x_u, U', các thách thức u_1, u_2,
    véc-tơ s, b_f = <s, b> = s(z) theo công thức tích, G_f = <s, G>, phép kiểm cuối.
(3) Bộ trích xuất của một vòng (n = 2, a = (17, 23), z = 5): ba nhánh thách thức
    2, 3, 5 -> hệ số beta -> khôi phục a.
(4) Phá binding trên nhóm toy: log rời rạc của G_i theo G1 -> mở C ra véc-tơ khác.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.field import GF
from snarklib.transcript import Transcript
from snarklib.pedersen import msm, smul, dlog_bruteforce
from snarklib.curve import toy_generator
from snarklib import ipa

F = GF(449)

# ---------------------------------------------------------------- (1)
n, z = 4, F(5)
pp = ipa.setup(n)
a = [F(c) for c in (3, 1, 4, 1)]
b = ipa.powers(z, n)
v = ipa.inner(a, b)
P = msm(a, pp.G) + v * pp.U
L, R, cL, cR = ipa.round_messages(pp.G, a, b, pp.U)
t = Transcript(b"ch4-mot-vong")
t.absorb(b"L", L); t.absorb(b"R", R)
x = ipa.challenge_nonzero(t, b"u")
xi = x.inverse()
print("(1) x =", x, " x^-1 =", xi, " x^2 =", x * x, " x^-2 =", (x * x).inverse())
print("    P = <a,G> + v.U =", P)
print("    a'_0 = x.3 + x^-1.4 =", x * 3 + xi * 4, "  a'_1 = x.1 + x^-1.1 =", x + xi)
print("    b'_0 = x^-1.1 + x.25 =", xi + x * 25, "  b'_1 = x^-1.5 + x.125 =", xi * 5 + x * 125)
G2, a2, b2 = ipa.fold_vectors(pp.G, a, b, x)
P2 = ipa.fold_commitment(P, L, R, x)
print("    P' = x^2 L + P + x^-2 R =", P2, "  <a',b'> =", ipa.inner(a2, b2),
      "  v + x^2 c_L + x^-2 c_R =", v + x * x * cL + (x * x).inverse() * cR)

# ---------------------------------------------------------------- (2)
C = ipa.commit(pp, a)
v1, proof = ipa.open(pp, a, z, Transcript(b"ch4-pcs"))
tv = Transcript(b"ch4-pcs")
tv.absorb(b"P", C); tv.absorb(b"z", z); tv.absorb(b"v", v1)
xu = ipa.challenge_nonzero(tv, b"x_u")
U1 = smul(pp.U, xu)
Q = C + v1 * U1
us = []
for Lj, Rj in zip(proof["L"], proof["R"]):
    tv.absorb(b"L", Lj); tv.absorb(b"R", Rj)
    u = ipa.challenge_nonzero(tv, b"u")
    us.append(u)
    Q = ipa.fold_commitment(Q, Lj, Rj, u)
s = ipa.s_vector(us)
bf = ipa.inner(s, b)
bf_prod = ipa.b_final_powers(us, z)
Gf = msm(s, pp.G)
af = F(proof["a"])
print("(2) x_u =", xu, " U' = x_u.U =", U1, " P' = C + v.U' =", C + v1 * U1)
print("    u =", us, " u^-1 =", [u.inverse() for u in us])
print("    s =", s, " b_f = <s,b> =", bf, " tích (u1^-1 + u1 z^2)(u2^-1 + u2 z) =", bf_prod)
print("    G_f = <s,G> =", Gf, " a =", af)
print("    Q =", Q, "  a.G_f + (a.b_f).U' =", af * Gf + (af * bf) * U1,
      "  khớp:", Q == af * Gf + (af * bf) * U1)
assert bf == bf_prod and Q == af * Gf + (af * bf) * U1

# ---------------------------------------------------------------- (3)
G = pp.G[:2]
U = pp.U
w = [F(17), F(23)]
bb = ipa.powers(z, 2)                       # (1, 5)
vv = ipa.inner(w, bb)
P0 = msm(w, G) + vv * U
L0, R0, _, _ = ipa.round_messages(G, w, bb, U)
xs = [F(2), F(3), F(5)]
leaves = []
for xc in xs:
    G1_, a1_, b1_ = ipa.fold_vectors(G, w, bb, xc)
    c = ipa.fold_commitment(P0, L0, R0, xc)
    assert c == a1_[0] * G1_[0] + (a1_[0] * b1_[0]) * U
    leaves.append((xc, a1_[0], b1_[0]))
print("(3) a = (17,23), b = (1,5), v =", vv, " P =", P0, " L =", L0, " R =", R0)
for xc, af_, bf_ in leaves:
    print("    nhánh x =", xc, ": a' =", af_, " b' =", bf_)
# giải hệ beta: sum beta = 1, sum beta x^2 = 0, sum beta x^-2 = 0
A = [[F(1), xc * xc, (xc * xc).inverse()] for xc in xs]
# beta là hàng đầu của A^{-1}: giải A^T beta = e_1
M = [[A[j][i] for j in range(3)] + [F(1) if i == 0 else F(0)] for i in range(3)]
for col in range(3):
    piv = next(r for r in range(col, 3) if M[r][col] != 0)
    M[col], M[piv] = M[piv], M[col]
    inv = M[col][col].inverse()
    M[col] = [e * inv for e in M[col]]
    for r in range(3):
        if r != col and M[r][col] != 0:
            f = M[r][col]
            M[r] = [er - f * ec for er, ec in zip(M[r], M[col])]
beta = [M[i][3] for i in range(3)]
print("    beta =", beta, " kiểm: sum =", sum(beta, F(0)),
      " sum beta x^2 =", sum((bt * xc * xc for bt, xc in zip(beta, xs)), F(0)),
      " sum beta x^-2 =", sum((bt * (xc * xc).inverse() for bt, xc in zip(beta, xs)), F(0)))
ext = [sum((bt * xc.inverse() * af_ for bt, (xc, af_, _) in zip(beta, leaves)), F(0)),
       sum((bt * xc * af_ for bt, (xc, af_, _) in zip(beta, leaves)), F(0))]
print("    véc-tơ trích xuất =", ext)
assert ext == w

# ---------------------------------------------------------------- (4)
_, _, G1 = toy_generator()
logs = [dlog_bruteforce(G1, Gi) for Gi in pp.G]
print("(4) log_{G1}(G_i) =", logs)
a_fake = [a[0] + logs[1], a[1] - logs[0], a[2], a[3]]
print("    a* = a + (g_1, -g_0, 0, 0) =", a_fake, " <a*,G> =", msm(a_fake, pp.G),
      " C =", C, " trùng:", msm(a_fake, pp.G) == C)
print("    p*(5) =", ipa.inner(a_fake, b), " khác p(5) =", v)
assert msm(a_fake, pp.G) == C
