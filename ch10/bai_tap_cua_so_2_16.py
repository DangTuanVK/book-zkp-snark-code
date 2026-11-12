"""Ch.10, đáp án bài tập: chọn cửa sổ Pippenger cho n = 2^16, b = 255 (đếm phép toán)."""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.pippenger import dem_pippenger, dem_ngay_tho, chi_phi_ly_thuyet, c_toi_uu_ly_thuyet

rng = random.Random(449)
n, b = 1 << 16, 255
s = [rng.getrandbits(b) for _ in range(n)]
print("n = 2^16, b = 255; c tối ưu theo công thức W(n + 2^(c+1)) + c(W-1): c =", c_toi_uu_ly_thuyet(n, b))
print("%3s %9s %8s %10s %10s" % ("c", "cộng", "nhân đôi", "tổng", "lý thuyết"))
tot = {}
for c in range(9, 16):
    d = dem_pippenger(s, c, b)
    tot[c] = d.tong()
    print("%3d %9d %8d %10d %10d" % (c, d.add, d.dbl, d.tong(), sum(chi_phi_ly_thuyet(n, b, c))))
c_best = min(tot, key=tot.get)
nt = dem_ngay_tho(s)
print("c tốt nhất (đếm thật) = %d; ngây thơ: %s, tổng %d -> gấp %.1f lần"
      % (c_best, nt, nt.tong(), nt.tong() / tot[c_best]))
assert c_best == 13
