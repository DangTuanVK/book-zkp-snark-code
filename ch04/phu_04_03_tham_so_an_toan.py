"""Chương 4 — Số truy vấn theo mức an toàn và sai số theo cỡ trường (tab:ch4-so-truy-van, tab:ch4-binius-tham-so).

(1) Số cột/truy vấn T nhỏ nhất với (1 - delta)^T <= 2^-lambda, tức
    T = ceil(lambda / log2(1 / (1 - delta))), cho bốn chế độ của delta:
      (a) delta -> (1 - rho)/3      (Ligero, Thaler Claim 10.6 / AHIV17 Lemma 4.2)
      (b) delta -> (1 - rho)/2      (bán kính giải mã duy nhất)
      (c) delta -> 1 - sqrt(rho)    (cận Johnson, Thaler Theorem 10.4 cho FRI, bỏ eta)
      (d) delta -> 1 - rho          (giả thuyết)
    Vì các cận là bất đẳng thức chặt, T ở giới hạn là cận dưới.
(2) Số hạng do trường của sum-check l*d/|F| (Chương 8) và |L0|/|F| của Ligero
    trên các trường nhị phân F_{2^k}, viết dưới dạng log2.
"""
import math

def T_of(lam, delta):
    return math.ceil(lam / math.log2(1 / (1 - delta)))

print("(1) Số truy vấn T theo rho, lambda")
print("%-6s %-26s %8s %8s %8s" % ("rho", "chế độ delta", "delta", "T(100)", "T(128)"))
rows = []
for rho in (1 / 2, 1 / 4, 1 / 16):
    modes = [("(a) (1-rho)/3", (1 - rho) / 3), ("(b) (1-rho)/2", (1 - rho) / 2),
             ("(c) 1-sqrt(rho)", 1 - math.sqrt(rho)), ("(d) 1-rho", 1 - rho)]
    for name, d in modes:
        t100, t128 = T_of(100, d), T_of(128, d)
        rows.append((rho, name, d, t100, t128))
        print("1/%-4d %-26s %8.4f %8d %8d" % (round(1 / rho), name, d, t100, t128))

# Kiểm chéo với Thaler: delta <= 1/3 cần hơn 200 lần lặp cho 2^-128
print("delta = 1/3, lambda = 128 -> T =", T_of(128, 1 / 3))
# Ví dụ số của Thaler (rho = 1/16, eta = 2^-14, delta = 3/4 - eta, 65 lần lặp)
eta = 2 ** -14
print("Thaler: (1 - (3/4 - eta))^65 = 2^%.2f" % (65 * math.log2(1 / 4 + eta)))

# Trên toy F_449, rho = 1/4, t = 4 và t = 8 (Mã nguồn 4.2 và Chương 8)
for t in (4, 8):
    print("toy: |L0|/|F| = 16/449 = %.4f ; (1 - 1/4)^%d = %.4f ; (1/4 + 1/4)^%d = %.4f"
          % (16 / 449, t, 0.75 ** t, t, 0.5 ** t))

print("(2) log2 của sai số do cỡ trường: sum-check l*d/|F| (l = 20, d = 3), Ligero |L0|/|F| (|L0| = 2^22)")
print("%-12s %14s %16s" % ("trường", "l*d/|F|", "|L0|/|F|"))
for k in (8, 16, 32, 64, 128):
    sc = math.log2(20 * 3) - k
    lg = 22 - k
    print("F_2^%-8d %14.2f %16.2f" % (k, sc, lg))

# (3) Kích thước proof Ligero (công thức (2k + T m) F + T log2(k/rho) H, 32 byte/phần tử)
#     với T theo chế độ (a) và (d) ở rho = 1/4, lambda = 100, cho n = 2^20 và 2^25.
print("(3) Proof Ligero, rho = 1/4, lambda = 100, 32 byte/phần tử")
for l in (20, 25):
    m, k = 1 << (l // 2), 1 << (l - l // 2)
    for name, T in (("(a)", T_of(100, 0.25)), ("(d)", T_of(100, 0.75))):
        nf, nh = 2 * k + T * m, T * int(math.log2(4 * k))
        print("n = 2^%d (m = %d, k = %d), T = %d %s: %d phần tử trường + %d băm = %d byte (%.2f MiB)"
              % (l, m, k, T, name, nf, nh, 32 * (nf + nh), 32 * (nf + nh) / 2 ** 20))
