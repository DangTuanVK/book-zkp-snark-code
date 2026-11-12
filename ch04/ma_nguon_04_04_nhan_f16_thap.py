"""Chương 4 — Nhân tay trong F_16 của tháp nhị phân (Mục 4.2d).

T_1 = F_4  = F_2[X0]/(X0^2 + X0 + 1)
T_2 = F_16 = F_4[X1]/(X1^2 + X0·X1 + 1)
Phần tử 4 bit b3 b2 b1 b0 = b0 + b1·X0 + b2·X1 + b3·X0·X1.
Ví dụ: a = 1 + X0 + X0·X1 (0b1011), b = X0 + X0·X1 (0b1010).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.gf2tower import mul, inv, fmt, mul_trace, T

# Bảng nhân F_4 (phần tử 0, 1, X0, X0+1 mã hóa 0..3)
print("Bảng nhân F_4:")
for x in range(4):
    print("  %-8s|" % fmt(x, 1), " ".join("%-8s" % fmt(mul(x, y, 1), 1) for y in range(4)))

a, b = 0b1011, 0b1010
tr = mul_trace(a, b, 2)                    # tách a = a0 + a1·X1, b = b0 + b1·X1
print("a =", fmt(a, 2), "  b =", fmt(b, 2))
print("a0 =", fmt(tr["a0"], 1), " a1 =", fmt(tr["a1"], 1),
      " b0 =", fmt(tr["b0"], 1), " b1 =", fmt(tr["b1"], 1), " (thuộc F_4)")
for key in ("a0b0", "a1b1", "a0b1+a1b0", "a1b1*c"):
    print("  %-10s = %s" % (key, fmt(tr[key], 1)))
print("Hệ số tự do = a0b0 + a1b1          =", fmt(tr["lo"], 1))
print("Hệ số của X1 = a0b1 + a1b0 + a1b1·X0 =", fmt(tr["hi"], 1))
ab = mul(a, b, 2)
print("a·b =", fmt(ab, 2), " (0b{:04b})".format(ab))

# Nghịch đảo của a và kiểm a·a^{-1} = 1
ai = inv(a, 2)
print("a^-1 =", fmt(ai, 2), " (0b{:04b});  a·a^-1 =".format(ai), fmt(mul(a, ai, 2), 2))

# Cùng chuỗi bit, nhân trong F_256 cho cùng kết quả (tháp lồng nhau)
F256 = T(3)
print("Trong F_256:", F256(a) * F256(b), " = a·b nhúng vào F_256")
assert ab == tr["ket_qua"] == 0b1110 and mul(a, ai, 2) == 1
assert (F256(a) * F256(b)).v == ab
