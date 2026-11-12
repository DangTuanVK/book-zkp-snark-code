"""Chương 4 — Một vòng gập đôi của IPA và PCS từ IPA (lst:ch4-ipa).

Đa thức p(X) = 3 + X + 4X^2 + X^3 trên F_449, mở tại z = 5.
Nhóm: nhóm con bậc 449 của đường cong toy (ký hiệu cộng, không dùng pairing).
Véc-tơ cơ sở G_0..G_3, U sinh tất định bằng băm-vào-đường-cong (snarklib.pedersen).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.field import GF
from snarklib.transcript import Transcript
from snarklib.pedersen import msm
from snarklib import ipa

F = GF(449)
n, z = 4, F(5)
pp = ipa.setup(n)                          # G_0..G_3, H, U "không ai biết log"
a = [F(c) for c in (3, 1, 4, 1)]           # hệ số của p(X)
b = ipa.powers(z, n)                       # b = (1, z, z^2, z^3)
v = ipa.inner(a, b)                        # v = p(z)
print("G =", pp.G, " U =", pp.U)
print("b =", b, " v = p(5) =", v)

# --- Một vòng gập đôi, viết tường minh ---
P = msm(a, pp.G) + v * pp.U                # P' = <a,G> + <a,b>·U
L, R, cL, cR = ipa.round_messages(pp.G, a, b, pp.U)
print("c_L = <a_lo,b_hi> =", cL, " c_R = <a_hi,b_lo> =", cR)
print("L =", L, " R =", R)
t = Transcript(b"ch4-mot-vong")
t.absorb(b"L", L); t.absorb(b"R", R)
x = ipa.challenge_nonzero(t, b"u")         # thách thức phụ thuộc L, R
G2, a2, b2 = ipa.fold_vectors(pp.G, a, b, x)
P2 = ipa.fold_commitment(P, L, R, x)       # x^2 L + P + x^-2 R
print("x =", x, " a' =", a2, " b' =", b2, " G' =", G2)
ok_round = P2 == msm(a2, G2) + ipa.inner(a2, b2) * pp.U
print("Đẳng thức sau một vòng:", ok_round)
assert ok_round

# --- PCS đầy đủ: log2 n = 2 vòng, Fiat–Shamir ---
C = ipa.commit(pp, a)
v1, proof = ipa.open(pp, a, z, Transcript(b"ch4-pcs"))
cnt = ipa.VCount()
ok = ipa.verify(pp, C, z, v1, proof, Transcript(b"ch4-pcs"), cnt=cnt)
print("C =", C, " proof: L =", proof["L"], " R =", proof["R"], " a =", proof["a"])
print("Kích thước proof:", ipa.proof_size_pcs(n), " Verifier:", cnt.as_dict())
print("Verifier chấp nhận:", ok)
bad = ipa.verify(pp, C, z, v1 + 1, proof, Transcript(b"ch4-pcs"))
print("Khai sai p(5) =", v1 + 1, "-> chấp nhận:", bad)
assert ok and v1 == v and not bad
