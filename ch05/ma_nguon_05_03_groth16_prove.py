"""Chương 5 -- Proof Groth16 thật cho N = 391 trên đường cong toy (Type 1).

Đường cong y^2 = x^3 + x trên F_5387, nhóm bậc r = 449, pairing đối xứng:
"G2" cũng là G1 = (2658, 1349). Seed cố định để output tái lập.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.curve import toy_generator
from snarklib.dsl import mach_phan_tich
from snarklib.groth16 import r1cs_to_qap, setup, prove, verify, combine

E, P, G = toy_generator()
g1 = g2 = G                                   # Type 1: G1 = G2
c = mach_phan_tich()[0]
r1cs = c.compile_r1cs()
qap = r1cs_to_qap(r1cs)                        # 13 + 2 hàng -> n = 16
x, w = c.r1cs_assignment({"N": 391, "a": 17, "b": 23})
z = [1] + x + w

rng = random.Random(2016)
pk, vk, tau = setup(qap, g1, g2, rng)
print("miền H: n = %d, ω = %d, t(X) = X^%d - 1" % (qap.n, int(qap.F.root_of_unity(qap.n)), qap.n))
print("tau (chỉ để minh họa) α, β, γ, δ, x =", [int(v) for v in tau])
print("CRS: |σ1| = %d điểm G1, |σ2| = %d điểm G2"
      % (3 + len(pk.xs1) + len(vk.IC1) + len(pk.K1) + len(pk.H1), 3 + len(pk.xs2)))
Ap, Bp, Cp, h, _ = combine(qap, z)
print("bậc A(X), B(X), C(X), h(X) =", Ap.degree(), Bp.degree(), Cp.degree(), h.degree())

proof = prove(pk, qap, z, rng)
print("A =", proof.A, " (G1)")
print("B =", proof.B, " (G2 = G1 trên toy)")
print("C =", proof.C, " (G1)")
ok = verify(vk, [391], proof)
sai = verify(vk, [392], proof)
print("verify(N = 391):", ok)
print("verify(N = 392):", sai)
assert ok and not sai
