"""Ch.10: nghịch đảo hàng loạt (mẹo Montgomery): n nghịch đảo trong F_p bằng
1 nghịch đảo + 3(n - 1) phép nhân. Ví dụ tay trên F_449 (n = 8) và đo trên p của BLS12-381."""
import os, sys, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ho_duong_cong import ARK_P_BLS12_381


def nghich_dao_hang_loat(a, p, dem):
    """Trả về [a_i^{-1} mod p]; mọi a_i phải khác 0. dem đếm 'nhan' và 'nghich_dao'."""
    n = len(a)
    tich = [a[0] % p]                         # tich[i] = a_0·a_1·...·a_i
    for i in range(1, n):
        tich.append(tich[-1] * a[i] % p); dem["nhan"] += 1
    t = pow(tich[-1], -1, p); dem["nghich_dao"] += 1   # t = (a_0···a_{n-1})^{-1}
    kq = [0] * n
    for i in range(n - 1, 0, -1):             # quét ngược
        kq[i] = tich[i - 1] * t % p            # a_i^{-1} = (a_0···a_{i-1})·t
        t = t * a[i] % p                       # t = (a_0···a_{i-1})^{-1}
        dem["nhan"] += 2
    kq[0] = t
    return kq


# Ví dụ tay: F_449, n = 8
p, a = 449, [2, 3, 5, 7, 11, 13, 17, 23]
dem = {"nhan": 0, "nghich_dao": 0}
inv = nghich_dao_hang_loat(a, p, dem)
print("F_449, a =", a)
print("a^{-1} =", inv)
print("chi phí: %d phép nhân + %d nghịch đảo (3(n-1) = %d)" % (dem["nhan"], dem["nghich_dao"], 3 * (len(a) - 1)))
assert all(x * y % p == 1 for x, y in zip(a, inv)) and dem["nhan"] == 3 * (len(a) - 1)

# Đo trên p của BLS12-381 (381 bit), n = 2^12, Python thuần
p = ARK_P_BLS12_381
rng = random.Random(449)
n = 1 << 12
a = [rng.randrange(1, p) for _ in range(n)]
t0 = time.perf_counter(); rieng = [pow(x, -1, p) for x in a]; t_rieng = time.perf_counter() - t0
dem = {"nhan": 0, "nghich_dao": 0}
t0 = time.perf_counter(); lo = nghich_dao_hang_loat(a, p, dem); t_lo = time.perf_counter() - t0
assert lo == rieng
t0 = time.perf_counter()
for i in range(n):
    _ = a[i] * a[i - 1] % p
t_nhan = (time.perf_counter() - t0) / n
t0 = time.perf_counter()
for i in range(256):
    _ = pow(a[i], -1, p)
t_inv = (time.perf_counter() - t0) / 256
print("\nBLS12-381 p (381 bit), n = %d:" % n)
print("  từng phần tử: %d nghịch đảo, %.1f ms" % (n, t_rieng * 1e3))
print("  hàng loạt   : %d nghịch đảo + %d phép nhân, %.1f ms (nhanh hơn %.1f lần)"
      % (dem["nghich_dao"], dem["nhan"], t_lo * 1e3, t_rieng / t_lo))
print("  một nghịch đảo ≈ %.0f phép nhân mod p (Python: %.2f µs / %.3f µs)"
      % (t_inv / t_nhan, t_inv * 1e6, t_nhan * 1e6))
print("  (thời gian phụ thuộc máy và tải; tỉ lệ mới là điều cần chú ý)")
