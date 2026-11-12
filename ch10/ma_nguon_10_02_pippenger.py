"""Ch.10 (lst:ch10-pippenger): MSM Pippenger (phương pháp xô, cửa sổ c) so với
nhân ngây thơ trên đường cong toy y^2 = x^3 + x / F_5387 (Type 1, G1 = G2; ở đây
chỉ dùng G1), kèm bộ đếm phép cộng / nhân đôi điểm.
Phần 1: kiểm đúng trên điểm thật. Phần 2: đếm phép toán với vô hướng b = 255 bit
(độ dài vô hướng của F_r BLS12-381) cho n = 2^4..2^12 -- số đếm không phụ thuộc đường cong."""
import os, sys, random
GOC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, GOC)
from snarklib.curve import toy_generator, R_TOY
from snarklib.pippenger import (NhomDiem, msm_pippenger, msm_ngay_tho, dem_pippenger,
                                dem_ngay_tho, chi_phi_ly_thuyet, c_toi_uu_ly_thuyet)

rng = random.Random(449)                       # tất định
E, P, G1 = toy_generator()

# ---- Phần 1: kiểm đúng, n = 16 điểm P_i = k_i·G1 ----
n = 16
k = [rng.randrange(1, R_TOY) for _ in range(n)]
diem = [ki * G1 for ki in k]
for b in (9, 255):                             # vô hướng trong F_449 và vô hướng 255 bit
    s = [rng.getrandbits(b) for _ in range(n)]
    mong_doi = (sum(si * ki for si, ki in zip(s, k)) % R_TOY) * G1
    Gp, Gn = NhomDiem(), NhomDiem()
    kq_p = msm_pippenger(Gp, diem, s, c=3, bits=b)
    kq_n = msm_ngay_tho(Gn, diem, s)
    assert kq_p == kq_n == mong_doi
    print("n = %d, b = %3d: MSM = %s | Pippenger c=3: %s | ngây thơ: %s"
          % (n, b, kq_p, Gp.dem, Gn.dem))

# ---- Phần 2: số phép toán theo n, b = 255, chọn c tốt nhất bằng thực nghiệm ----
b = 255
print("\nb = %d bit; 'tổng' = cộng + nhân đôi; LT = W(n + 2^(c+1)) + c(W-1)" % b)
print("%6s %3s %3s %9s %9s %9s %10s %8s" % ("n", "c*", "cLT", "Pip.cộng", "Pip.đôi",
                                         "LT(c*)", "ngây thơ", "tỉ lệ"))
dong_csv = ["n,c_tot_nhat,c_ly_thuyet,pip_cong,pip_nhan_doi,ly_thuyet,ngay_tho,ti_le"]
for logn in range(4, 13):
    n = 1 << logn
    s = [rng.getrandbits(b) for _ in range(n)]
    dem = {c: dem_pippenger(s, c, b) for c in range(1, 15)}
    c_best = min(dem, key=lambda c: dem[c].tong())
    d, nt = dem[c_best], dem_ngay_tho(s)
    lt = sum(chi_phi_ly_thuyet(n, b, c_best))
    ti_le = nt.tong() / d.tong()
    print("%6d %3d %3d %9d %9d %9d %10d %7.1fx" % (n, c_best, c_toi_uu_ly_thuyet(n, b),
                                                  d.add, d.dbl, lt, nt.tong(), ti_le))
    dong_csv.append("%d,%d,%d,%d,%d,%d,%d,%.2f" % (n, c_best, c_toi_uu_ly_thuyet(n, b),
                                                  d.add, d.dbl, lt, nt.tong(), ti_le))
    assert d.tong() <= lt and ti_le > 1

# Độ nhạy theo c ở n = 2^12
print("\nn = 4096: tổng phép toán theo c =",
      ", ".join("%d:%d" % (c, dem[c].tong()) for c in range(4, 13)))
os.makedirs(os.path.join(GOC, "ket_qua", "ch10"), exist_ok=True)
with open(os.path.join(GOC, "ket_qua", "ch10", "pippenger.csv"), "w") as f:
    f.write("\n".join(dong_csv) + "\n")
