"""Chương 4 — Kiểm các tính chất của tháp nhị phân dùng trong Mục 4.2.

(1) Đa thức X^2 + X_{i-1}·X + 1 không có nghiệm trong T_i (i = 0..4), nên
    T_{i+1} = T_i[X_i]/(X_i^2 + X_{i-1}X_i + 1) là một trường (đa thức bậc 2
    không nghiệm thì bất khả quy). Kiểm bằng vét cạn.
(2) Mọi phần tử khác 0 của T_i thỏa a^(|T_i| - 1) = 1 (i = 1..3, vét cạn).
(3) Phần tử của F_2 = T_0 thỏa x^2 = x: ràng buộc bit b(b - 1) = 0 tự động đúng;
    trong T_3 chỉ 0 và 1 thỏa x^2 = x.
(4) Nhân "nhỏ với lớn": phần tử T_1 nhân phần tử T_3 bằng 4 phép nhân trong T_1
    (T_3 là không gian véc-tơ chiều 4 trên T_1, cơ sở 1, X1, X2, X1·X2).
(5) Đóng gói: 2^20 bit witness, mỗi bit một phần tử 32 byte vs 128 bit/phần tử T_7.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.gf2tower import mul, gen, power, bits_of

for lvl in range(0, 5):
    c = gen(lvl)                       # X_{lvl-1} (bằng 1 khi lvl = 0)
    size = 1 << bits_of(lvl)
    roots = [x for x in range(size) if mul(x, x, lvl) ^ mul(c, x, lvl) ^ 1 == 0]
    print("T_%d (%6d phần tử): số nghiệm của X^2 + X_%d·X + 1 = %d"
          % (lvl, size, lvl - 1, len(roots)))
    assert not roots

for lvl in (1, 2, 3):
    q = 1 << bits_of(lvl)
    assert all(power(a, q - 1, lvl) == 1 for a in range(1, q))
    print("T_%d: a^%d = 1 với mọi a khác 0 (vét cạn %d phần tử)" % (lvl, q - 1, q - 1))

idem3 = [x for x in range(256) if mul(x, x, 3) == x]
print("Trong F_2: 0^2 = 0, 1^2 = 1;  trong T_3 các x với x^2 = x:", idem3)
assert idem3 == [0, 1]

# Nhân T_1 x T_3 theo từng tọa độ (4 khối 2 bit) so với nhân đầy đủ trong T_3
s, big = 0b10, 0b10110110            # s = X0 thuộc T_1
by_coord = 0
for i in range(4):
    by_coord |= mul(s, (big >> (2 * i)) & 3, 1) << (2 * i)
print("X0 · 0b10110110: nhân trong T_3 = 0b{:08b}, theo 4 tọa độ T_1 = 0b{:08b}"
      .format(mul(s, big, 3), by_coord))
assert by_coord == mul(s, big, 3)

n_bits = 1 << 20
print("2^20 bit: mỗi bit 32 byte -> %d byte; đóng gói 128 bit/phần tử T_7 -> %d phần tử, %d byte"
      % (n_bits * 32, n_bits // 128, (n_bits // 128) * 16))
