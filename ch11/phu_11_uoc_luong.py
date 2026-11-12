"""Ch.11 (mã phụ, không in sách): ngân sách chu kỳ cho MSM và NTT.

Phần 1: xung đột bucket -- mô phỏng một bộ cộng điểm có pipeline sâu L, mỗi chu
        kỳ phát một phép cộng "bucket[k] += P_i" theo thứ tự; nếu bucket k còn một
        phép cộng chưa ra khỏi pipeline thì phải chờ (dừng in-order) hoặc gác lại
        (hàng chờ). Đếm số chu kỳ trên số phép cộng, so với xác suất lý thuyết.
Phần 2: ước lượng chu kỳ và thời gian cho MSM 2^16..2^20 (vô hướng 254 bit) theo
        eq:ch11-chu-ky-msm và NTT theo eq:ch11-chu-ky-ntt, dùng tần số/diện tích
        DC lsi_10k của ket_qua/ch11/tong_hop.csv (montmul_w256 tổng hợp với r của
        BLS12-381, dùng làm đại diện cho p 254 bit của BN254). ĐÂY LÀ ƯỚC LƯỢNG.
"""
import os, sys, random, csv, math
GOC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, GOC)

rng = random.Random(11)

# ---- Phần 1: xung đột bucket ----
def mo_phong(c, L, so_phep, kieu):
    """Trả về số chu kỳ để phát hết so_phep phép cộng với 2^c - 1 bucket."""
    B = (1 << c) - 1
    xong = {}                               # bucket -> chu kỳ kết quả ra khỏi pipeline
    chu_ky, cho = 0, []
    nguon = (rng.randrange(1, B + 1) for _ in range(so_phep))
    con = so_phep
    while con > 0:
        k = None
        if kieu == "hang_cho":              # ưu tiên phép đang gác nếu bucket đã rảnh
            for idx, kk in enumerate(cho):
                if xong.get(kk, -1) <= chu_ky:
                    k = cho.pop(idx)
                    break
            if k is None and len(cho) < 64:
                kk = next(nguon, None)
                if kk is not None:
                    if xong.get(kk, -1) <= chu_ky:
                        k = kk
                    else:
                        cho.append(kk)
        else:                               # dừng in-order
            if not cho:
                cho.append(next(nguon))
            if xong.get(cho[0], -1) <= chu_ky:
                k = cho.pop(0)
        if k is not None:
            xong[k] = chu_ky + L
            con -= 1
        chu_ky += 1
    return chu_ky

print("Phần 1: xung đột bucket, 2^16 phép cộng ngẫu nhiên đều trên 2^c - 1 bucket")
print("%3s %4s %12s %12s %14s" % ("c", "L", "dừng in-order", "hàng chờ 64", "1-(1-1/B)^(L-1)"))
for c in (8, 12, 16):
    for L in (8, 32):
        B = (1 << c) - 1
        n = 1 << 16
        d = mo_phong(c, L, n, "dung")
        h = mo_phong(c, L, n, "hang_cho")
        p = 1 - (1 - 1 / B) ** (L - 1)
        print("%3d %4d %12.3f %12.3f %14.3f" % (c, L, d / n, h / n, p))

# ---- Phần 2: ước lượng ----
dc = {}
with open(os.path.join(GOC, "ket_qua", "ch11", "tong_hop.csv")) as f:
    for d in csv.DictReader(f):
        if d["chon"] == "1":
            dc[d["top"]] = (float(d["f_uoc_luong_MHz"]), float(d["dien_tich"]))

b = 254                                  # vô hướng BN254 (r 254 bit)
print("\nPhần 2a: MSM, b = %d bit, số phép cộng A(n,c) = ceil(b/c)·(n + 2^(c+1)), "
      "chu kỳ ~ A/P" % b)
print("%6s %3s %14s %16s %16s" % ("n", "c", "A (cộng điểm)", "P=1 chu kỳ", "P=16 chu kỳ"))
for logn in (16, 18, 20):
    n = 1 << logn
    best = min(range(4, 21), key=lambda c: math.ceil(b / c) * (n + (1 << (c + 1))))
    A = math.ceil(b / best) * (n + (1 << (best + 1)))
    print("%6s %3d %14d %16d %16d" % ("2^%d" % logn, best, A, A, math.ceil(A / 16)))

# một đơn vị cộng điểm XYZZ madd-2008-s: 8M + 2S = 10 phép nhân trường (EFD)
f256, a256 = dc["montmul_w256"]
n = 1 << 20
best = min(range(4, 21), key=lambda c: math.ceil(b / c) * (n + (1 << (c + 1))))
A = math.ceil(b / best) * (n + (1 << (best + 1)))
print("\nPhần 2b: đơn vị MSM giả định cho n = 2^20 (BN254 G1: p 254 bit vừa W = 256; lsi_10k, f >= %.2f MHz)" % f256)
print("  diện tích 1 montmul W=256 = %.0f; bộ cộng điểm = 10 montmul (pipeline, 1 phép cộng/chu kỳ)"
      % a256)
for P in (1, 4, 16):
    cyc = math.ceil(A / P)
    print("  P = %2d bộ cộng: %d montmul, diện tích ~ %.2e, chu kỳ ~ %d, thời gian ~ %.2f s ở %.2f MHz"
          % (P, 10 * P, 10 * P * a256, cyc, cyc / (f256 * 1e6), f256))

print("\nPhần 2c: NTT n = 2^20, (n/2)·log2 n butterfly, B butterfly song song")
fbf, abf = dc["ntt_bfly_gl64"]
n = 1 << 20
nb = (n // 2) * 20
for B in (1, 4, 16):
    cyc = math.ceil(nb / B)
    print("  Goldilocks W=64, B = %2d: %d butterfly, chu kỳ ~ %d, thời gian ~ %.1f ms ở %.2f MHz, diện tích ~ %.2e"
          % (B, nb, cyc, 1e3 * cyc / (fbf * 1e6), fbf, B * abf))
