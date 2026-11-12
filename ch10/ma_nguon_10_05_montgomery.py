"""Ch.10: nhân Montgomery (REDC). Ví dụ tay trên p = 449, R = 2^9 với 17·23;
vét cạn mọi cặp trong F_449; kiểm trên p của BLS12-381 với R = 2^384 (6 từ 64 bit, như arkworks)."""
import os, sys, random
GOC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, GOC); sys.path.insert(0, os.path.join(GOC, "ch10"))
from snarklib.montgomery import Montgomery
from ho_duong_cong import ARK_P_BLS12_381

# 1) Ví dụ tay: p = 449, R = 2^9 = 512
M = Montgomery(449)
print("p = %d, R = 2^%d = %d, p' = -p^{-1} mod R = %d, R mod p = %d, R^2 mod p = %d"
      % (M.p, M.k, M.R, M.p_prime, M.R_mod_p, M.R2))
a, b = 17, 23
am, bm = M.to_mont(a), M.to_mont(b)
T = am * bm
m = ((T % M.R) * M.p_prime) % M.R
t = (T + m * M.p) // M.R
print("ã = 17·R mod p = %d, b̃ = 23·R mod p = %d" % (am, bm))
print("T = ã·b̃ = %d; m = (T mod R)·p' mod R = %d; t = (T + m·p)/R = %d%s"
      % (T, m, t, " -> t - p = %d" % (t - M.p) if t >= M.p else " (< p, không trừ)"))
cm = M.mul(am, bm)
print("REDC(T) = %d = 391·R mod p; from_mont -> %d" % (cm, M.from_mont(cm)))
assert cm == 391 * M.R % 449 and M.from_mont(cm) == 391

# 2) Vét cạn F_449: mọi cặp (a, b), 449^2 = 201601 phép nhân
M.so_lan_tru_cuoi = 0
loi = sum(1 for x in range(449) for y in range(449)
          if M.from_mont(M.mul(M.to_mont(x), M.to_mont(y))) != x * y % 449)
print("vét cạn 449^2 = %d cặp: %d lỗi" % (449 * 449, loi))
assert loi == 0

# 3) BLS12-381: p 381 bit, R = 2^384 (s = 6 từ 64 bit), CIOS đếm phép nhân từ máy
M = Montgomery(ARK_P_BLS12_381, word=64)
rng = random.Random(449)
M.so_lan_tru_cuoi = 0
for _ in range(2000):
    x, y = rng.randrange(M.p), rng.randrange(M.p)
    xm, ym = M.to_mont(x), M.to_mont(y)
    z, nhan = M.mul_cios(xm, ym)
    assert z == M.mul(xm, ym) and M.from_mont(z) == x * y % M.p
print("\n%s, s = %d từ: 2000 phép nhân ngẫu nhiên khớp a·b mod p" % (M, M.s))
print("CIOS: %d phép nhân từ máy 64x64 (2s^2 + s = %d)" % (nhan, 2 * M.s**2 + M.s))
print("tỉ lệ REDC cần bước trừ cuối: %.1f%% (trên 8000 lần REDC: 2 to_mont + mul + from_mont mỗi vòng)"
      % (100 * M.so_lan_tru_cuoi / 8000))
assert nhan == 2 * M.s**2 + M.s
