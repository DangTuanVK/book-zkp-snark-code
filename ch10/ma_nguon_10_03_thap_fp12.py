"""Ch.10: tháp F_p -> F_{p^2} -> F_{p^6} -> F_{p^12} trên số nguyên tố nhỏ p = 19.
Chọn p = 19 vì là số nguyên tố nhỏ nhất có cùng đồng dư với p của BLS12-381
(p ≡ 3 mod 4 để u^2 + 1 bất khả quy; p ≡ 1 mod 3, tức p ≡ 7 mod 12) và với
ξ = 1 + u (như BLS12-381) vừa không chính phương vừa không lập phương trong F_{p^2}."""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thap_truong import Thap, DEM, xi_hop_le

# 1) Dò p ≡ 7 (mod 12) nhỏ nhất dùng được với ξ = 1 + u
for p in (7, 19):
    khong_cp, khong_lp = xi_hop_le(p)
    print("p = %2d: ξ không chính phương = %s, không lập phương = %s" % (p, khong_cp, khong_lp))
p = 19
T = Thap(p)
assert xi_hop_le(p) == (True, True) and xi_hop_le(7)[0] is False

# 2) Đếm phép nhân F_p cho một phép nhân F_{p^12}
rng = random.Random(449)
a, b = T.ngau_nhien12(rng), T.ngau_nhien12(rng)
DEM["mul"] = 0; ab = T.f12_mul(a, b);                  n_kara = DEM["mul"]
DEM["mul"] = 0; ab2 = T.f12_mul(a, b, karatsuba=False); n_truong = DEM["mul"]
assert ab == ab2
print("nhân F_{p^12}: Karatsuba %d phép nhân F_p, trường học %d" % (n_kara, n_truong))
DEM["mul"] = DEM["inv"] = 0; ai = T.f12_inv(a)
print("nghịch đảo F_{p^12}: %d phép nhân F_p + %d nghịch đảo F_p" % (DEM["mul"], DEM["inv"]))

# 3) Kiểm tiên đề trường trên 200 bộ (x, y, z) ngẫu nhiên tất định
one, zero = T.one12(), T.zero12()
M, A = T.f12_mul, T.f12_add
for _ in range(200):
    x, y, z = (T.ngau_nhien12(rng) for _ in range(3))
    assert M(x, y) == M(y, x) and M(M(x, y), z) == M(x, M(y, z))
    assert M(x, A(y, z)) == A(M(x, y), M(x, z)) and M(x, one) == x
    if x != zero:
        assert M(x, T.f12_inv(x)) == one
print("tiên đề (giao hoán, kết hợp, phân phối, đơn vị, nghịch đảo): đạt trên 200 mẫu")

# 4) Cấu trúc: w^2 = v, w^6 = ξ; Frobenius x -> x^p; w có bậc 12 trên F_p
w = T.w()
w6 = T.f12_pow(w, 6)
print("w^6 =", T.he_so(w6), "(= ξ = 1 + u)")
assert w6 == (((1, 1), (0, 0), (0, 0)), ((0, 0),) * 3)
k = next(k for k in range(1, 13) if T.f12_pow(w, p**k) == w)
print("k nhỏ nhất với w^(p^k) = w: k = %d -> [F_p(w) : F_p] = 12" % k)
x = T.ngau_nhien12(rng)
assert T.f12_pow(x, p**12) == x and T.f12_pow(x, p**12 - 1) == one and k == 12
assert T.f12_pow(M(a, b), p) == M(T.f12_pow(a, p), T.f12_pow(b, p))
print("x^(p^12) = x, x^(p^12 - 1) = 1, (ab)^p = a^p b^p: đạt; |F_{p^12}| = 19^12 = %d" % p**12)
