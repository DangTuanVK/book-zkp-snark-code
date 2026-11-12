"""Chương 5 -- Biên dịch mạch ví dụ N = 391 bằng DSL tí hon sang R1CS.

Mệnh đề: biết a, b với a·b = N, a và b vừa 5 bit. Điều kiện a, b >= 2 không
cần cổng riêng: a = 0 cho a·b = 0, a = 1 buộc b = 391 > 31 (xem 05_02).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.dsl import Circuit, mach_phan_tich

c = Circuit()                       # trường F_449
N = c.public("N")                   # đầu vào công khai
a = c.private("a")                  # witness
b = c.private("b")
c.to_bits(a, 5, prefix="a")         # 5 ràng buộc bit + 1 ràng buộc tổ hợp
c.to_bits(b, 5, prefix="b")
c.assert_mul(a, b, N)               # a·b = N

r1cs = c.compile_r1cs()
order = [c.names[v] for v in c.z_order()]
print("z =", "(" + ", ".join(order) + ")")
print("R1CS: m = %d ràng buộc, n = %d biến (1 hằng + %d công khai + %d witness)"
      % (r1cs.m, r1cs.n, r1cs.n_public, r1cs.n - 1 - r1cs.n_public))
for j, (A, B, C, nhan) in enumerate(c.constraints):
    print("  %2d %-9s (%s) * (%s) = (%s)" % (j, nhan, A, B, C))

x, w = c.r1cs_assignment({"N": 391, "a": 17, "b": 23})
print("x =", [int(v) for v in x])
print("w =", [int(v) for v in w])
print("thỏa R1CS (17, 23):", r1cs.is_satisfied(x, w))
x2, w2 = c.r1cs_assignment({"N": 391, "a": 17, "b": 24})
print("thỏa R1CS (17, 24):", r1cs.is_satisfied(x2, w2))

# Tự kiểm: trùng với mạch dựng sẵn trong snarklib.dsl
ref = mach_phan_tich()[0].compile_r1cs()
assert (r1cs.A, r1cs.B, r1cs.C) == (ref.A, ref.B, ref.C)
assert r1cs.is_satisfied(x, w) and not r1cs.is_satisfied(x2, w2)
assert (r1cs.m, r1cs.n) == (13, 14)
