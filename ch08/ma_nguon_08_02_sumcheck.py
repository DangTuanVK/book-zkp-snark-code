"""Chương 8 — Sum-check cho tích hai MLE (lst:ch8-sumcheck).

H = Σ_{b in {0,1}^2} f~(b)·h~(b),  f = (17, 23, 28, 30),  h = (23, 17, 30, 28)
  = 17·23 + 23·17 + 28·30 + 30·28 = 2·391 + 2·840 ≡ 4·391 ≡ 217 (mod 449).
(840 ≡ 391 là nghiệm giả do tràn trường của ví dụ xuyên suốt.)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.field import GF
from snarklib.transcript import Transcript
from snarklib import sumcheck as sc, mle

F = GF(449)
f, h = [17, 23, 28, 30], [23, 17, 30, 28]
g = lambda u: u[0] * u[1]                  # g(u, v) = u·v, bậc 2 theo mỗi biến

proof, r, finals = sc.prove(F, [f, h], g, 2, Transcript(b"ch8-sumcheck"))
print("Khẳng định H =", proof["claim"])
claim = proof["claim"]
for j, s in enumerate(proof["polys"]):
    # Prover gửi s_j(0), s_j(1), s_j(2); Verifier kiểm s_j(0) + s_j(1) = khẳng định hiện tại
    print("Vòng %d: s_%d(0..2) = %s, s(0)+s(1) = %d (cần %d), r_%d = %d, s_%d(r_%d) = %d"
          % (j + 1, j + 1, s, s[0] + s[1], claim, j + 1, r[j], j + 1, j + 1,
             sc.interpolate_eval(F, s, r[j])))
    claim = sc.interpolate_eval(F, s, r[j])
print("Bảng sau khi gập x1 = r1:", mle.fold(F, f, r[0]), mle.fold(F, h, r[0]))
print("Oracle: f~(r) =", finals[0], ", h~(r) =", finals[1], ", tích =", g(finals))

ok, rr, final = sc.verify(F, 2, 2, proof, Transcript(b"ch8-sumcheck"))
oracle = [mle.evaluate(F, f, rr), mle.evaluate(F, h, rr)]
print("Verifier: kiểm vòng", ok, "; kiểm cuối g(oracle(r)) =", g(oracle), "= s_2(r_2) =", final)
print("Proof: %d vòng x %d giá trị = %d phần tử trường; cận lỗi l·d/|F| = 4/449 = %.4f"
      % (2, 3, 6, sc.soundness_bound(2, 2, 449)))

# Khẳng định sai H' = 218 với đa thức trung thực: bị bác ngay ở vòng 1
bad = dict(proof, claim=F(218))
print("Khẳng định sai 218 ->", sc.verify(F, 2, 2, bad, Transcript(b"ch8-sumcheck"))[0])
assert proof["claim"] == 217 == (2 * 391 + 2 * 840) % 449
assert ok and g(oracle) == final and not sc.verify(F, 2, 2, bad, Transcript(b"ch8-sumcheck"))[0]
