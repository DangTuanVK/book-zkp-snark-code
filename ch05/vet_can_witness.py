"""Vét cạn (a, b) ∈ [0, 31]^2 thỏa MẠCH ví dụ N = 391 (dành cho Chương 12).

Bit và tổ hợp buộc a, b ∈ [0, 31], nên vét cạn trên [0, 31]^2 là đủ: mọi
witness thỏa đều có (a, b) trong miền này. Kiểm cả ba bản biên dịch và Groth16.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.curve import toy_generator
from snarklib.dsl import mach_phan_tich
from snarklib.groth16 import r1cs_to_qap, setup, prove, verify

c = mach_phan_tich()[0]
r1cs, pk_c, ccs = c.compile_r1cs(), c.compile_plonkish(), c.compile_ccs()
nghiem = []
for a in range(32):
    for b in range(32):
        inp = {"N": 391, "a": a, "b": b}
        x, w = c.r1cs_assignment(inp)
        ok = [r1cs.is_satisfied(x, w), pk_c.is_satisfied(*c.plonkish_assignment(inp)),
              ccs.is_satisfied(x, w)]
        assert len(set(ok)) == 1
        if ok[0]:
            nghiem.append((a, b, a * b))
print("Số cặp đã thử: 1024; số witness thỏa mạch:", len(nghiem))
for a, b, p in nghiem:
    print("  (a, b) = (%2d, %2d): a·b = %3d = %d·449 + %d" % (a, b, p, p // 449, p % 449))

E, P, G = toy_generator()
qap = r1cs_to_qap(r1cs)
pk, vk, _ = setup(qap, G, G, random.Random(2016))
for a, b, _ in nghiem:
    x, w = c.r1cs_assignment({"N": 391, "a": a, "b": b})
    print("  Groth16 proof cho (%d, %d): verify =" % (a, b),
          verify(vk, x, prove(pk, qap, [1] + x + w, random.Random(a * 32 + b))))
assert [(a, b) for a, b, _ in nghiem] == [(17, 23), (23, 17), (28, 30), (30, 28)]
