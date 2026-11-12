"""Ch.11: mô hình bit-chính-xác của montmul.sv (W = 16, R = 2^16, p = 449).
In giá trị từng dây T, m, U, t, y cho ví dụ 17·23 -- đúng tên tín hiệu trong RTL --
và đếm số lần bước trừ cuối xảy ra khi vét cạn 449^2 cặp (so với R = 2^9 ở Ch.10)."""
import os, sys
GOC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, GOC)
from snarklib.montgomery import Montgomery

W, P = 16, 449
M = Montgomery(P, word=W)                 # R = 2^16 (một từ 16 bit)
R, PINV = M.R, M.p_prime


def montmul_rtl(a, b):
    """Từng dòng ứng với một dòng always_comb của montmul.sv."""
    T = a * b                             # T = a_q * b_q        (2W bit)
    m = (T & (R - 1)) * PINV & (R - 1)    # m = T[W-1:0] * PINV  (cắt W bit)
    U = T + m * P                         # U = {1'b0,T} + m*P   (2W+1 bit)
    t = U >> W                            # t = U[2W:W]
    y = t - P if t >= P else t            # trừ có điều kiện
    return T, m, U, t, y


print("W = %d, R = 2^%d = %d, p = %d" % (W, W, R, P))
print("PINV = -p^{-1} mod R = %d = 0x%04x;  R mod p = %d;  R^2 mod p = %d"
      % (PINV, PINV, M.R_mod_p, M.R2))
a, b = M.to_mont(17), M.to_mont(23)       # dạng Montgomery: 17R mod p, 23R mod p
T, m, U, t, y = montmul_rtl(a, b)
print("a = 17R mod p = %d (0x%04x),  b = 23R mod p = %d (0x%04x)" % (a, a, b, b))
print("T = %d,  m = %d,  U = %d = %d·2^16,  t = %d,  y = %d" % (T, m, U, U >> W, t, y))
print("U mod R = %d (chia hết);  y = 391·R mod p = %d;  from_mont(y) = %d"
      % (U % R, 391 * R % P, M.from_mont(y)))
assert U % R == 0 and y == M.mul(a, b) == 391 * R % P and M.from_mont(y) == 391
print("vector RTL (a b y, hex): %04x%04x%04x" % (a, b, y))

# Vét cạn: bước trừ cuối xảy ra bao nhiêu lần? (cũng là số lỗi khi bỏ bước trừ trong RTL)
for k in (9, 16):
    Mk = Montgomery(P, word=k)
    dem = sum(1 for x in range(P) for z in range(P)
              if ((x * z + ((x * z * Mk.p_prime) & Mk.mask) * P) >> k) >= P)
    print("R = 2^%-2d: bước trừ cuối %6d / %d cặp = %.3f%%" % (k, dem, P * P, 100 * dem / P / P))
    if k == 16:
        assert dem == 155                  # khớp số lỗi của bản RTL bỏ bước trừ (VCS)
loi = sum(1 for x in range(P) for z in range(P) if montmul_rtl(x, z)[4] != M.mul(x, z))
print("mô hình RTL vs snarklib.montgomery: %d cặp, %d lỗi" % (P * P, loi))
assert loi == 0
