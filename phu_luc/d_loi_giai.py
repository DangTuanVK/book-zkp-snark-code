"""Phụ lục D -- Kiểm mọi đáp số của phần lời giải bài tập chọn lọc bằng snarklib.

Mỗi khối tương ứng một bài tập; các phép `assert` chốt đúng những giá trị được in
trong Phụ lục D. Khối Bài tập 12.3 cần z3 (gói z3-solver).
Chạy từ thư mục gốc kho mã:  python phu_luc/d_loi_giai.py
"""
import math
import os
import random
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from snarklib import GF, gf2tower as T, kzg
from snarklib import uniqueness as U
from snarklib.curve import Curve, toy_generator
from snarklib.dsl import Circuit, mach_phan_tich, mach_phan_tich_sua
from snarklib.groth16 import prove, r1cs_to_qap, setup, verify
from snarklib.montgomery import Montgomery
from snarklib.poly import Poly, domain, lagrange_interpolate, vanishing_poly

F = GF(449)
P_BN254 = 21888242871839275222246405745257275088696311157297823662689037894645226208583


def tieu_de(s):
    print("\n=== %s ===" % s)


# ---------------------------------------------------------------- Bài tập 1.2
tieu_de("Bài tập 1.2: omega_32")
g = F(3)
day = [int(g ** e) for e in (1, 2, 4, 8)]
w32 = g ** 14
print("3^1, 3^2, 3^4, 3^8 =", day)
print("3^8 * 3^4 =", int(g ** 12), "; 3^12 * 3^2 =", int(w32))
print("omega_32^16 =", int(w32 ** 16), "(= -1); omega_32^2 =", int(w32 ** 2))
assert day == [3, 9, 81, 275] and int(g ** 12) == 274 and int(w32) == 221
assert int(w32 ** 16) == 448 and int(w32 ** 2) == 349 and w32.order() == 32
H32, H16 = set(map(int, domain(F, 32))), set(map(int, domain(F, 16)))
assert H16 < H32
print("H_16 nằm trong H_32:", H16 < H32)

# ---------------------------------------------------------------- Bài tập 1.4
tieu_de("Bài tập 1.4: y^2 = x^3 + x trên F_2693")
E2693 = Curve(2693, 1, 0)
so_diem = E2693.count_points()
print("q mod 4 =", 2693 % 4, "; #E =", so_diem, "; q + 1 =", 2694,
      "; vết t =", 2694 - so_diem)
print("449 | #E:", so_diem % 449 == 0, "; #E mod 449 =", so_diem % 449)
print("phân tích #E:", {p: k for p, k in __import__("snarklib").factorize(so_diem).items()})
assert 2693 % 4 == 1 and so_diem != 2694
assert so_diem == 2788 == 2 ** 2 * 17 * 41 and so_diem % 449 == 94

# ---------------------------------------------------------------- Bài tập 2.3
tieu_de("Bài tập 2.3: khoảng cách thống kê của mặt nạ rút gọn modulo 449")


def kc_rut_gon(bits, p=449):
    s = 2 ** bits % p
    return s, Fraction(s * (p - s), 2 ** bits * p)


s16, d16 = kc_rut_gon(16)
kiem = Fraction(1, 2) * (s16 * abs(Fraction(2 ** 16 // 449 + 1, 2 ** 16) - Fraction(1, 449))
                         + (449 - s16) * abs(Fraction(2 ** 16 // 449, 2 ** 16) - Fraction(1, 449)))
assert d16 == kiem
print("2^16 = %d*449 + %d; Delta = %s = %.4e; log2(Delta) = %.2f" % (
    2 ** 16 // 449, s16, d16, float(d16), math.log2(d16)))
s512, d512 = kc_rut_gon(512)
print("2^512 mod 449 = %d; log2(Delta) = %.2f" % (
    s512, math.log2(d512.numerator) - math.log2(d512.denominator)))
assert (s16, 2 ** 16 // 449) == (431, 145)
assert d16 == Fraction(3879, 14712832) and s512 == 359
assert round(math.log2(d16), 2) == -11.89
assert round(math.log2(d512.numerator) - math.log2(d512.denominator), 2) == -505.83

# ---------------------------------------------------------------- Bài tập 3.1
tieu_de("Bài tập 3.1: gamma công bố trước khi gửi y_i")
Fk = kzg.F
rng = random.Random(3)
srs = kzg.setup(15, 77)
fs = [Poly(Fk, [rng.randrange(449) for _ in range(16)]) for _ in range(4)]
Cs = [kzg.commit(srs, f) for f in fs]
z, gam = Fk(7), Fk(5)
ys = [f(z) for f in fs]
fg = Poly.zero(Fk)
for i, f in enumerate(fs):
    fg = fg + f * (gam ** i)
v, W = kzg.open(srs, fg, z)
Fgop = Cs[0]
for i in range(1, 4):
    Fgop = Fgop + int(gam ** i) * Cs[i]


def kiem_gop(ds):
    vv = sum((gam ** i * y for i, y in enumerate(ds)), Fk(0))
    return kzg.verify(srs, Fgop, z, vv, W)


ys_gia = list(ys)
ys_gia[0] = ys[0] + 1
ys_gia[1] = ys[1] - gam.inverse()
print("y =", [int(y) for y in ys], "; y' =", [int(y) for y in ys_gia], "; W =", W)
print("gamma^-1 =", int(gam.inverse()), "; trung thực:", kiem_gop(ys), "; khẳng định sai:", kiem_gop(ys_gia))
ys_chi_sua_1 = [ys[0] + 1] + ys[1:]
print("chỉ sửa y_1 (không bù y_2):", kiem_gop(ys_chi_sua_1))
assert kiem_gop(ys) and kiem_gop(ys_gia) and not kiem_gop(ys_chi_sua_1)
assert ys_gia[0] != fs[0](z) and ys_gia[1] != fs[1](z)
assert [int(y) for y in ys] == [202, 58, 433, 292] and [int(y) for y in ys_gia] == [203, 417, 433, 292]
assert str(W) == "(381, 4532)" and int(gam.inverse()) == 90

# ---------------------------------------------------------------- Bài tập 4.3
tieu_de("Bài tập 4.3: nhân và nghịch đảo trong F_256")
a43 = 0b00010011              # 1 + X0 + X2
b43 = 0b10000010              # X0 + X0 X1 X2
vet = T.mul_trace(a43, b43, 3)
print("a =", T.fmt(a43, 3), "=", format(a43, "08b"))
print("b =", T.fmt(b43, 3), "=", format(b43, "08b"))
for k in ("a0", "a1", "b0", "b1", "c", "a0b0", "a1b1", "a0b1+a1b0", "a1b1*c", "lo", "hi"):
    print("  %-10s = %-22s (%s)" % (k, T.fmt(vet[k], 2), format(vet[k], "04b")))
ab = T.mul(a43, b43, 3)
ainv = T.inv(a43, 3)
print("ab =", T.fmt(ab, 3), "=", format(ab, "08b"))
print("a^-1 =", T.fmt(ainv, 3), "=", format(ainv, "08b"), "; a*a^-1 =", T.mul(a43, ainv, 3))
x2_vuong = T.mul(0x10, 0x10, 3)
x2_ba = T.mul(x2_vuong, 0x10, 3)
print("X2^2 =", T.fmt(x2_vuong, 3), "; X2^3 =", T.fmt(x2_ba, 3))
assert vet["ket_qua"] == ab and T.mul(a43, ainv, 3) == 1
a0, a1, c_t = vet["a0"], vet["a1"], vet["c"]
delta = T.mul(a0, a0, 2) ^ T.mul(T.mul(a0, a1, 2), c_t, 2) ^ T.mul(a1, a1, 2)
print("Delta =", T.fmt(delta, 2), "; Delta^-1 =", T.fmt(T.inv(delta, 2), 2))
assert delta == 0b1111 and T.inv(delta, 2) == 0b0111
assert x2_vuong == T.add(1, T.mul(0b100, 0x10, 3))
assert ab == 0b10001001 and ainv == 0b01111011 and x2_ba == 0b10000100
_goc = T.mul
dem_and = [0]


def _mul_dem(a, b, lvl):
    if lvl == 0:
        dem_and[0] += 1
    return _goc(a, b, lvl)


T.mul = _mul_dem
T.mul(a43, b43, 3)
T.mul = _goc
print("số phép AND của mul ở T_3:", dem_and[0])
assert dem_and[0] == 64


def mul_c_tuyen_tinh(a, lvl):
    """Nhân với hằng c = X_{lvl-1} như ánh xạ F_2-tuyến tính: chỉ XOR, không AND."""
    ket_qua = 0
    for i in range(1 << lvl):
        if a >> i & 1:
            ket_qua ^= _goc(1 << i, T.gen(lvl), lvl)
    return ket_qua


def mul_27(a, b, lvl, dem):
    if lvl == 0:
        dem[0] += 1
        return a & b
    h = 1 << (lvl - 1)
    mask = (1 << h) - 1
    a0, a1, b0, b1 = a & mask, a >> h, b & mask, b >> h
    p00, p11 = mul_27(a0, b0, lvl - 1, dem), mul_27(a1, b1, lvl - 1, dem)
    pm = mul_27(a0 ^ a1, b0 ^ b1, lvl - 1, dem) ^ p00 ^ p11
    return (p00 ^ p11) | ((pm ^ mul_c_tuyen_tinh(p11, lvl - 1)) << h)


dem27 = [0]
assert all(mul_27(u, v, 3, dem27) == _goc(u, v, 3) for u in range(256) for v in range(256))
print("nhân với c bằng XOR: %d phép AND mỗi phép nhân, đúng trên 256^2 cặp" % (dem27[0] // 256 ** 2))
assert dem27[0] == 27 * 256 ** 2

# ---------------------------------------------------------------- Bài tập 5.4
tieu_de("Bài tập 5.4: đổi ví dụ sang N = 437")
c437 = mach_phan_tich()[0]
r437 = c437.compile_r1cs()
pl437 = c437.compile_plonkish()
print("R1CS m = %d, Plonkish %d hàng (không đổi so với N = 391)" % (r437.m, pl437.n))
ng437 = []
for a in range(32):
    for b in range(32):
        x, w = c437.r1cs_assignment({"N": 437, "a": a, "b": b})
        if r437.is_satisfied(x, w):
            ng437.append((a, b))
print("witness trong [0,31]^2:", ng437, "; k trong ab = 437 + 449k:",
      sorted({(a * b - 437) // 449 for a, b in ng437}))
print("437 + 449 = 886 = 2*443; ước <= 31 của 886:", [d for d in range(1, 32) if 886 % d == 0])
assert ng437 == [(19, 23), (23, 19)]
assert r437.m == 13 and pl437.n == 20
E, _, G = toy_generator()
qap = r1cs_to_qap(r437)
rng = random.Random(2016)
pk, vk, tau = setup(qap, G, G, rng)
for a, b in ng437:
    x, w = c437.r1cs_assignment({"N": 437, "a": a, "b": b})
    pi = prove(pk, qap, [1] + x + w, random.Random(a * 32 + b))
    ok = verify(vk, x, pi)
    print("  Groth16 proof cho (%d, %d): verify = %s" % (a, b, ok))
    assert ok
x, w = c437.r1cs_assignment({"N": 437, "a": 19, "b": 23})
pi = prove(pk, qap, [1] + x + w, rng)
al, be, ga, de, xx = tau
L1 = be * qap.u[1](xx) + al * qap.v[1](xx) + qap.w[1](xx)
K = int(L1 / de) * G
dem_gia = sum(verify(vk, [Nm], pi._replace(C=pi.C + int(qap.F(437 - Nm)) * K)) for Nm in range(449))
print("số N' có proof giả từ proof của N = 437:", dem_gia)
assert dem_gia == 449

# ---------------------------------------------------------------- Bài tập 6.2
tieu_de("Bài tập 6.2: cổng tùy biến c = a^5")
n = 32
Hn = domain(F, n)
ZH = vanishing_poly(F, n)
rng = random.Random(62)
giatri_a = [F(rng.randrange(449)) for _ in Hn]
aH = lagrange_interpolate(F, Hn, giatri_a)
cH = lagrange_interpolate(F, Hn, [v ** 5 for v in giatri_a])
qP = lagrange_interpolate(F, Hn, [F(rng.randrange(1, 449)) for _ in Hn])
aX = aH + Poly(F, [F(rng.randrange(1, 449)), F(rng.randrange(1, 449))]) * ZH
cX = cH + Poly(F, [F(rng.randrange(1, 449)), F(rng.randrange(1, 449))]) * ZH
tu = qP * (aX ** 5 - cX)
t62, du = tu.divmod(ZH)
print("deg q_P = %d, deg a = %d, deg tử số = %d (6n+4 = %d), deg t = %d (5n+4 = %d), dư = 0: %s"
      % (qP.degree(), aX.degree(), tu.degree(), 6 * n + 4, t62.degree(), 5 * n + 4, du.is_zero()))
assert tu.degree() == 196 and t62.degree() == 164 and du.is_zero()
so_phan = (t62.degree() - (n + 4)) // n + 1
print("số phần (phần cuối bậc n + 4 = %d): %d; thêm %d phần tử G1 so với 3" % (n + 4, so_phan, so_phan - 3))
assert so_phan == 5
so_coset = -(-(t62.degree() + 1) // n)
sinh = F(3)
dai_dien = [sinh ** j for j in range(1, 14)]    # 13 coset khác H (F*/H cyclic cấp 14)
assert all(k ** n != 1 for k in dai_dien)
assert len({int(k ** n) for k in dai_dien}) == 13
diem = [k * h for k in dai_dien[:so_coset] for h in Hn]
gia_tri = [tu(u) / ZH(u) for u in diem]
assert lagrange_interpolate(F, diem[:165], gia_tri[:165]) == t62
print("cần %d điểm -> %d coset cỡ %d; số coset khác H của F*/H: %d; nội suy lại t: đúng"
      % (t62.degree() + 1, so_coset, n, 448 // n - 1))
assert so_coset == 6 and 448 // n - 1 == 13

# ---------------------------------------------------------------- Bài tập 7.2
tieu_de("Bài tập 7.2: Lasso cho bảng 16 bit")


def khuc(v, c=4, l=4):
    return [(v >> (l * (c - 1 - i))) & (2 ** l - 1) for i in range(c)]


def g4(y):
    return sum(2 ** (4 * (3 - i)) * u for i, u in enumerate(y))


def mle_bang_con(x):
    return sum(2 ** (4 - s) * u for s, u in enumerate(x, 1))


for v in (17, 23, 391, 840):
    kh = khuc(v)
    assert g4(kh) == v
    bit = [[(u >> (3 - j)) & 1 for j in range(4)] for u in kh]
    assert [mle_bang_con(b) for b in bit] == kh
    print("  %3d -> khúc %s, g = %d" % (v, kh, g4(kh)))
nl, d = 1000, 2 ** 16
chi_phi = {4: 3 * 4 * nl + 4 * 2 ** 4, 2: 3 * 2 * nl + 2 * 2 ** 8, 1: nl + d, 16: 3 * 16 * nl + 16 * 2}
print("số phần tử cam kết (c: giá trị):", chi_phi)
assert chi_phi == {4: 12064, 2: 6512, 1: 66536, 16: 48032}
print("phân tích bit (7.bit) cho 1000 giá trị 16 bit: R1CS %d ràng buộc, %d biến bit" % (17 * nl, 16 * nl))

# ---------------------------------------------------------------- Bài tập 8.2
tieu_de("Bài tập 8.2: sum-check bậc hai")
f_mle = lambda x1, x2: F(17) + 11 * x1 + 6 * x2 - 4 * x1 * x2
g82 = lambda x1, x2: f_mle(x1, x2) ** 2
B = [F(0), F(1)]
H82 = sum((g82(u, v) for u in B for v in B), F(0))
s1 = lambda X: sum((g82(X, v) for v in B), F(0))
kiem_s1 = lambda X: (F(17) + 11 * X) ** 2 + (F(23) + 7 * X) ** 2
assert all(s1(F(u)) == kiem_s1(F(u)) for u in range(449))
r1, r2 = F(2), F(3)
s2 = lambda X: g82(r1, X)
print("H = %d, s1(0) = %d, s1(1) = %d, s1(2) = %d" % (int(H82), int(s1(F(0))), int(s1(F(1))), int(s1(r1))))
print("f~(2,X) = 39 - 2X: f~(2,0) = %d, f~(2,1) = %d; s2(0) = %d, s2(1) = %d, s2(3) = %d"
      % (int(f_mle(r1, F(0))), int(f_mle(r1, F(1))), int(s2(F(0))), int(s2(F(1))), int(s2(r2))))
assert (int(H82), int(s1(F(0))), int(s1(F(1)))) == (257, 369, 337)
assert s2(F(0)) + s2(F(1)) == s1(r1) and int(s2(r2)) == 191 and int(f_mle(r1, r2)) == 33
assert (int(s1(r1)), int(s2(F(0))), int(s2(F(1)))) == (196, 174, 22) and int(kiem_s1(F(1))) == 337
kD = F(8).inverse()
D = Poly(F, [F(6), F(-5), F(1)]) * kD        # k(X - 2)(X - 3)
assert D(F(0)) + D(F(1)) == 1
print("k = 8^-1 = %d; D(X) = k(X-2)(X-3); D(0) + D(1) = %d" % (int(kD), int(D(F(0)) + D(F(1)))))


def thoat(ri1, ri2):
    """Chiến lược hai vòng: lệch còn lại ở vòng 2 được đặt vào D2 có hai nghiệm 2, 3."""
    lech = D(ri1)
    if lech == 0:
        return True
    D2 = Poly(F, [F(6), F(-5), F(1)]) * (lech * F(8).inverse())
    assert D2(F(0)) + D2(F(1)) == lech
    return D2(ri2) == 0


dem_thoat = sum(thoat(F(u), F(v)) for u in range(449) for v in range(449))
print("số cặp (r1, r2) thoát: %d / %d = %.4f%%; cận 4/449 = %.4f%%; chỉ riêng vòng 1: 2/449 = %.4f%%"
      % (dem_thoat, 449 ** 2, 100 * dem_thoat / 449 ** 2, 400 / 449, 200 / 449))
assert dem_thoat == 1792

# ---------------------------------------------------------------- Bài tập 9.3
tieu_de("Bài tập 9.3: CCS bậc 3 và số hạng chéo")
f93 = lambda u, y, x: x ** 3 - u ** 2 * y
z1, z2 = (1, 8, 2), (1, 27, 3)
u1, y1, x1 = z1
u2, y2, x2 = z2
e1 = 3 * x1 ** 2 * x2 - (2 * u1 * u2 * y1 + u1 ** 2 * y2)
e2 = 3 * x1 * x2 ** 2 - (u2 ** 2 * y1 + 2 * u1 * u2 * y2)
zg = tuple(p + 5 * q for p, q in zip(z1, z2))
print("f(z1) = %d, f(z2) = %d, e1 = %d, e2 = %d, z1 + 5 z2 = %s, f = %d, mod 449 = %d, 5e1 + 25e2 = %d"
      % (f93(*z1), f93(*z2), e1, e2, zg, f93(*zg), f93(*zg) % 449, 5 * e1 + 25 * e2))
assert (e1, e2, zg, f93(*zg)) == (-7, -8, (6, 143, 17), -235)
assert f93(*zg) % 449 == 214 == (5 * e1 + 25 * e2) % 449
for rr in range(449):
    zr = tuple(p + rr * q for p, q in zip(z1, z2))
    assert (f93(*zr) - (rr * e1 + rr ** 2 * e2)) % 449 == 0
print("khai triển đúng với mọi r trong F_449")

# ---------------------------------------------------------------- Bài tập 10.3
tieu_de("Bài tập 10.3: Montgomery cho 28 * 30")
M = Montgomery(449)
assert M.R == 512
at, bt = 28 * 512 % 449, 30 * 512 % 449
Tm = at * bt
m = ((Tm & 511) * M.p_prime) & 511
t = (Tm + m * 449) >> 9
print("a~ = %d, b~ = %d, T = %d, m = %d, t = %d, REDC(t) = %d" % (at, bt, Tm, m, t, M.redc(t)))
assert (at, bt, Tm, m, t) == (417, 94, 39198, 354, 387) and M.mul(at, bt) == 387
assert M.mul(17 * 512 % 449, 23 * 512 % 449) == 387 and M.redc(387) == 391
Mb = Montgomery(P_BN254, word=64)
am, bm = Mb.to_mont(28), Mb.to_mont(30)
kq, so_nhan = Mb.mul_cios(am, bm)
print("BN254: s = %d, số phép nhân từ của CIOS = %d (2s^2 + s = %d), 28*30 = %d"
      % (Mb.s, so_nhan, 2 * Mb.s ** 2 + Mb.s, Mb.from_mont(kq)))
assert so_nhan == 36 and Mb.from_mont(kq) == 840

# ---------------------------------------------------------------- Bài tập 11.4
tieu_de("Bài tập 11.4: rút gọn lười")
p = 449
R16 = 1 << 16
pp = (-pow(p, -1, R16)) % R16


def redc_luoi(Tv):
    mv = ((Tv & (R16 - 1)) * pp) & (R16 - 1)
    return (Tv + mv * p) >> 16


max_t = 0
for u in range(2 * p):
    for v in range(2 * p):
        tv = redc_luoi(u * v)
        assert (tv * R16 - u * v) % p == 0
        max_t = max(max_t, tv)
print("R = 2^16, 4p = %d, R > 4p: %s; max t trên [0,2p)^2 = %d < 2p = %d" % (4 * p, R16 > 4 * p, max_t, 2 * p))
assert R16 > 4 * p and max_t < 2 * p
assert max_t == 460
loi = [(u, v) for u in range(p) for v in range(p) if redc_luoi(u * v) >= p]
gia_tri_t = sorted({redc_luoi(u * v) for u, v in loi})
print("số cặp [0,p)^2 có t >= p: %d; các giá trị t: %s" % (len(loi), gia_tri_t))
assert len(loi) == 155 and all((redc_luoi(u * v) - p) == (u * v * pow(R16, -1, p)) % p for u, v in loi)
assert gia_tri_t == [450, 451]
rngc = random.Random(114)
for _ in range(2000):
    x_luoi = rngc.randrange(p)
    x_dung = x_luoi
    kq_luoi = (x_luoi * R16) % p
    kq_dung = kq_luoi
    for _ in range(20):
        y = rngc.randrange(p)
        y_m = y * R16 % p
        kq_luoi = redc_luoi(kq_luoi * y_m)
        assert kq_luoi < 2 * p
        x_dung = x_dung * y % p
    assert redc_luoi(kq_luoi) % p == x_dung
print("2000 chuỗi 20 phép nhân không trừ: luôn < 2p, rút gọn một lần ở cuối cho đúng tích")

# ---------------------------------------------------------------- Bài tập 12.2
tieu_de("Bài tập 12.2: so sánh bằng bit, không dùng gadget")


def mach_bit():
    c = Circuit()
    vN, a, b = c.public("N"), c.private("a"), c.private("b")
    ba = c.to_bits(a, 5, prefix="a")
    bb = c.to_bits(b, 5, prefix="b")
    c.enforce(ba[4], ba[3], c.const(0), "le_bit")
    c.enforce(ba[4], ba[2], c.const(0), "le_bit")
    c.enforce(bb[4], bb[3], c.const(0), "le_bit")
    c.assert_mul(a, b, vN)
    return c


for x in range(32):
    a4, a3, a2 = (x >> 4) & 1, (x >> 3) & 1, (x >> 2) & 1
    assert (x <= 19) == (a4 * a3 == 0 and a4 * a2 == 0)
    assert (x <= 23) == (a4 * a3 == 0)
cb = mach_bit()
rb, plb = cb.compile_r1cs(), cb.compile_plonkish()
ws = U.enumerate_witnesses(rb, [391])
print("mạch mới: R1CS m = %d, Plonkish %d hàng; witness tại N = 391: %s" % (
    rb.m, plb.n, [tuple(w[:2]) for w in ws]))
assert rb.m == 16 and [tuple(w[:2]) for w in ws] == [(17, 23)]
assert plb.n == 23
gia = 0
tong_w = 0
for Nv in range(1, 449):
    for w in U.enumerate_witnesses(rb, [Nv]):
        tong_w += 1
        ai, bi = int(w[0]), int(w[1])
        if ai * bi != Nv or ai > 19 or bi > 23:
            gia += 1
print("quét N trong [1, 448]: %d witness, %d nghiệm giả" % (tong_w, gia))
assert gia == 0
assert tong_w == 437
sua = mach_phan_tich_sua()[0]
print("so với gadget: R1CS %d, Plonkish %d hàng" % (sua.compile_r1cs().m, sua.compile_plonkish().n))
assert sua.compile_r1cs().m == 25 and sua.compile_plonkish().n == 38
mien = lambda c: (1 << (c.compile_plonkish().n - 1).bit_length(), r1cs_to_qap(c.compile_r1cs()).n)
print("miền (PLONK, QAP): mạch bit %s, gadget %s" % (mien(cb), mien(sua)))
assert mien(cb) == (32, 32) and mien(sua) == (64, 32)

# ---------------------------------------------------------------- Bài tập 12.3
tieu_de("Bài tập 12.3: mạch cho N = 437")
cu = mach_phan_tich()[0].compile_r1cs()
wc = [tuple(w[:2]) for w in U.enumerate_witnesses(cu, [437])]
print("W_C(437) của mạch 5 bit:", wc)
assert wc == [(19, 23), (23, 19)]
Ba = math.isqrt(437)
Bb = 448 // Ba
csua = mach_phan_tich_sua(bounds=(Ba, Bb))[0].compile_r1cs()
wsua = U.enumerate_witnesses(csua, [437])
print("cận theo công thức: B_a = %d, B_b = %d, B_a*B_b = %d; số witness = %d" % (Ba, Bb, Ba * Bb, len(wsua)))
assert (Ba, Bb) == (20, 22) and wsua == []
c2 = mach_phan_tich_sua(bounds=(19, 23))[0]
r2c = c2.compile_r1cs()
w2 = U.enumerate_witnesses(r2c, [437])
print("cận (19, 23): 19*23 = %d < 449; witness = %s" % (19 * 23, [tuple(w[:2]) for w in w2]))
assert [tuple(w[:2]) for w in w2] == [(19, 23)]
print("z3, witness thứ hai của cận (19, 23):", U.second_witness_z3(r2c, [437], w2[0]))
assert U.second_witness_z3(r2c, [437], w2[0]) is None

print("\nTất cả phép kiểm của Phụ lục D đều đạt.")
