"""Chương 8 — Đo xác suất lọt của Prover gian trong sum-check (thm:ch8-sumcheck).

Prover gian khẳng định H' = H + 1 (sai). Ở mỗi vòng j, Prover tính s_j trung thực
và gửi g_j = s_j + D_j với D_j(X) = δ_j·(X - 5)(X - 7)/59, δ_j = khẳng định hiện
tại - (s_j(0) + s_j(1)); nhờ đó g_j(0) + g_j(1) luôn khớp. Prover chỉ thoát khi
r_j ∈ {5, 7} (D_j(r_j) = 0, khẳng định trở về đúng). Xác suất lý thuyết
1 - (1 - 2/449)^2 ≈ 0,00889 ≤ l·d/|F| = 4/449.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.field import GF
from snarklib.transcript import Transcript
from snarklib import sumcheck as sc, mle

F = GF(449)
f, h = [17, 23, 28, 30], [23, 17, 30, 28]
g = lambda u: u[0] * u[1]
D = lambda x: (F(x) - 5) * (F(x) - 7) / 59      # D(0) + D(1) = (35 + 24)/59 = 1


def cheat(label):
    """Chạy song song Prover gian và Verifier trên một transcript; trả về Verifier có chấp nhận."""
    tp, tv = Transcript(label), Transcript(label)
    T = [[F(v) for v in f], [F(v) for v in h]]
    H = sum((g([T[0][i], T[1][i]]) for i in range(4)), F(0))
    claim = H + 1                                   # khẳng định sai
    polys, rs = [], []
    sc._absorb_header(F, tp, b"sc", 2, 2, claim)
    for j in range(2):
        half = len(T[0]) // 2
        s = [sum((g([t[i] + x * (t[i + half] - t[i]) for t in T]) for i in range(half)), F(0))
             for x in range(3)]
        delta = claim - (s[0] + s[1])
        gj = [s[x] + delta * D(x) for x in range(3)]
        polys.append(gj)
        tp.absorb(b"sc/s%d" % j, [int(v) for v in gj])
        rj = sc._challenge(F, tp, b"sc/r%d" % j)
        rs.append(rj)
        claim = sc.interpolate_eval(F, gj, rj)
        T = [[t[i] + rj * (t[i + half] - t[i]) for i in range(half)] for t in T]
    proof = {"claim": H + 1, "polys": polys}
    ok, r, final = sc.verify(F, 2, 2, proof, tv)
    return ok and g([mle.evaluate(F, f, r), mle.evaluate(F, h, r)]) == final


trials = 4000
acc = sum(cheat(b"gian-%d" % i) for i in range(trials))
theory = 1 - (1 - 2 / 449) ** 2
print("Số lần thử: %d; Verifier chấp nhận khẳng định sai: %d lần (tỉ lệ %.5f)" % (trials, acc, acc / trials))
print("Lý thuyết với chiến lược này: %.5f (kỳ vọng %.1f lần); cận l·d/|F| = %.5f"
      % (theory, theory * trials, sc.soundness_bound(2, 2, 449)))
assert 0 < acc < 3 * theory * trials
