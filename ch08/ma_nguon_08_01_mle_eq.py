"""Chương 8 — MLE và đa thức eq trên ví dụ nhỏ (Mục 8.1).

Bảng f = (17, 23, 28, 30) trên {0,1}^2 (big-endian: f(0,0), f(0,1), f(1,0), f(1,1)),
tức hai cặp thừa số (17, 23) và (28, 30) của ví dụ xuyên suốt.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.field import GF
from snarklib import mle
from snarklib.sumcheck import counting_field

F = GF(449)
f = [17, 23, 28, 30]

# Khai triển f~(x1, x2) = c0 + c1·x1 + c2·x2 + c12·x1·x2 (giải từ 4 giá trị)
c0 = F(f[0]); c2 = F(f[1]) - f[0]; c1 = F(f[2]) - f[0]; c12 = F(f[3]) - f[2] - f[1] + f[0]
print("f~(x1,x2) = %d + %d·x1 + %d·x2 + %d·x1·x2   (%d = -4 mod 449)" % (c0, c1, c2, c12, c12))

r = [F(2), F(3)]
print("eq_table(r = (2,3)) =", mle.eq_table(F, r))          # [eq(b, r) cho b = 00, 01, 10, 11]
v = mle.evaluate(F, f, r)
v_fold = mle.evaluate_by_folding(F, f, r)
v_coef = c0 + c1 * r[0] + c2 * r[1] + c12 * r[0] * r[1]
print("f~(2,3) theo eq =", v, "; theo gập đôi =", v_fold, "; theo hệ số =", v_coef)
print("Gập x1 = 2:", mle.fold(F, f, 2), "-> gập x2 = 3:", mle.fold(F, mle.fold(F, f, 2), 3))

# Một mở rộng KHÔNG đa tuyến tính: g = f~ + x1·(x1 - 1) trùng f trên {0,1}^2, khác tại r
g = lambda x: c0 + c1 * x[0] + c2 * x[1] + c12 * x[0] * x[1] + x[0] * (x[0] - 1)
print("g trên siêu lập phương:", [g([F(b1), F(b2)]) for b1 in (0, 1) for b2 in (0, 1)],
      "; g(2,3) =", g(r), "≠", v)

# Đánh giá O(2^l): đếm phép nhân cho l = 2..10 (eq_table + tích vô hướng)
print("l | số phép nhân của mle.evaluate | 3·2^l - 2")
CF = counting_field(F)
for l in range(2, 11, 2):
    CF.counter["mul"] = 0
    mle.evaluate(CF, [CF(i) for i in range(1 << l)], [CF(5 + i) for i in range(l)])
    print("%2d | %5d | %5d" % (l, CF.counter["mul"], 3 * (1 << l) - 2))
    assert CF.counter["mul"] == 3 * (1 << l) - 2

assert v == v_fold == v_coef and g(r) != v
assert [g([F(b1), F(b2)]) for b1 in (0, 1) for b2 in (0, 1)] == f
