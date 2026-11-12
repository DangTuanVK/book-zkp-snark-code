"""Chương 8 — GKR trên mạch 3 tầng tính tay được (Mục 8.2b).

Vào (a, b, c, d) = (17, 23, 28, 30) — hai cặp thừa số (17, 23) và (28, 30) của N = 391.
Tầng 1: [a·b, c·d, a+b, c+d];  tầng 0 (ra): [a·b + c·d, (a+b)·(c+d)].
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.transcript import Transcript
from snarklib import gkr, mle

circ, x = gkr.vi_du_mach()
F = circ.F
W = circ.evaluate(x)
for i, Wi in enumerate(W):
    print("W_%d =" % i, Wi, " (k_%d = %d)" % (i, circ.k[i]))
for i in range(circ.d):
    for typ in ("add", "mul"):
        T = circ.predicate_table(i, typ)
        ones = [mle.bits(j, circ.k[i] + 2 * circ.k[i + 1]) for j, v in enumerate(T) if v]
        print("%s_%d = 1 tại (a|b|c) =" % (typ, i), ones)

proof = gkr.prove(circ, x, Transcript(b"ch8-gkr"))
for i, L in enumerate(proof["layers"]):
    print("Tầng %d: khẳng định %d, sum-check %d biến (bậc 2), z1 = W~_%d(b*) = %d, z2 = W~_%d(c*) = %d, q =" %
          (i, L["claim"], 2 * circ.k[i + 1], i + 1, L["z"][0], i + 1, L["z"][1]), L["q"])
trace = []
ok = gkr.verify(circ, x, W[0], proof, Transcript(b"ch8-gkr"), trace=trace)
for i, r, v in trace:
    print("  khẳng định W~_%d(%s) = %d" % (i, r, v))
print("Verifier chấp nhận:", ok, "; proof:", gkr.proof_size(circ), "phần tử trường")

# Prover khai sai đầu ra (334 thay 333): thử 4 nhãn transcript, cho phép lọt tối đa 1 lần (~1/449)
bad = sum(gkr.verify(circ, x, [334, 75], proof, Transcript(b"gkr-sai-%d" % k)) for k in range(4))
print("Đầu ra sai [334, 75] được chấp nhận:", bad, "/ 4 lần")
assert ok and W[0] == [333, 75] and W[1] == [391, 391, 40, 58] and bad <= 1
