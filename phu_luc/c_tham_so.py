"""Phụ lục C (lst:appc-kiem): tự kiểm bảng tham số đường cong của sách.

Nguồn giá trị (chép nguyên văn, không tính lại bằng tay):
  BN254, BLS12-381: u, p, r từ arkworks 0.6.0 (ch10/ho_duong_cong.py);
                    BLS12-381 còn đối chiếu README của zkcrypto/bls12_381.
  BLS12-377:        u từ El Housni--Guillevic 2020, công thức (6);
                    p, r từ Bảng 2 của bài đó và Hình 16 của Zexe.
  BW6-761:          r, q, đa thức, h_t = 13, h_y = 9 từ Bảng 3, 4 của El Housni--Guillevic.
  Pallas, Vesta:    p, q từ README của zcash/pasta và zcash/pasta_curves.
Chạy: python phu_luc/c_tham_so.py  (đặt CARGO_HOME để đối chiếu thêm với crate).
"""
import glob, math, os, re, sys

GOC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [GOC, os.path.join(GOC, "ch10")]
from ho_duong_cong import bn, bls12, miller_rabin, two_adicity, U_BN254, U_BLS12_381
from snarklib.hang_so import compute

U_377 = 0x8508c00000000001
R_377 = 0x12ab655e9a2ca55660b44d1e5c37b00159aa76fed00000010a11800000000001
P_377 = int("1ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1"
            "ef3622fba094800170b5d44300000008508c00000000001", 16)
Q_BW6 = int("122e824fb83ce0ad187c94004faff3eb926186a81d14688528"
            "275ef8087be41707ba638e584e91903cebaff25b423048689c8e"
            "d12f9fd9071dcd3dc73ebff2e98a116c25667a8f8160cf8aeeaf"
            "0a437e6913e6870000082f49d00000000008b", 16)
P_PASTA = 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001
Q_PASTA = 0x40000000000000000000000000000000224698fc0994a8dd8c46eb2100000001


def tu_64_bit(n):
    """Viết hệ 16, tách thành từng từ 64 bit tính từ bit thấp."""
    s = "%x" % n
    dau = len(s) % 16
    khoi = ([s[:dau]] if dau else []) + [s[i:i + 16] for i in range(dau, len(s), 16)]
    return " ".join(khoi)


def can_bac_hai(a, p):
    """Tonelli--Shanks; trả về None nếu a không chính phương."""
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    s, m = two_adicity(p - 1), p - 1
    m >>= s
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    c, t, x, e = pow(z, m, p), pow(a, m, p), pow(a, (m + 1) // 2, p), s
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2, i = t2 * t2 % p, i + 1
        b = pow(c, 1 << (e - i - 1), p)
        x, c, t, e = x * b % p, b * b % p, t * b * b % p, i
    return x


def cong(P, Q, a, p):
    """Cộng affine trên y^2 = x^3 + ax + b; None là điểm vô cực."""
    if P is None or Q is None:
        return P or Q
    if P[0] == Q[0] and (P[1] + Q[1]) % p == 0:
        return None
    if P == Q:
        m = (3 * P[0] * P[0] + a) * pow(2 * P[1], -1, p) % p
    else:
        m = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, p) % p
    x = (m * m - P[0] - Q[0]) % p
    return x, (m * (P[0] - x) - P[1]) % p


def nhan(k, P, a, p):
    R = None
    while k:
        if k & 1:
            R = cong(R, P, a, p)
        P, k = cong(P, P, a, p), k >> 1
    return R


def cac_diem(b, p, so=3):
    """Ba điểm tất định: x nhỏ nhất có x^3 + b chính phương."""
    ds, x = [], 1
    while len(ds) < so:
        y = can_bac_hai(x**3 + b, p)
        if y:
            ds.append((x, y))
        x += 1
    return ds


def bac_nhom_j0(p, b, t):
    """Chứng minh #E = p + 1 - t cho E: y^2 = x^3 + b (j = 0, p ≡ 1 mod 3).

    Từ 4p = t^2 + 3y^2, #E chỉ có thể là một trong sáu giá trị p + 1 - t_i.
    Nhân ba điểm với sáu ứng viên; ứng viên duy nhất triệt tiêu cả ba là #E.
    """
    y2, du = divmod(4 * p - t * t, 3)
    y = math.isqrt(y2)
    assert du == 0 and y * y == y2 and p % 3 == 1
    v1, v2 = (t + 3 * y) // 2, (t - 3 * y) // 2
    vet = {t, -t, v1, v2, -v1, -v2}
    assert len(vet) == 6
    diem = cac_diem(b, p)
    con = [v for v in vet if all(nhan(p + 1 - v, P, 0, p) is None for P in diem)]
    assert len(con) == 1, "không xác định duy nhất được #E"
    return p + 1 - con[0]


def bac_nhung(p, r, gioi_han=10000):
    x = 1
    for k in range(1, gioi_han + 1):
        x = x * p % r
        if x == 1:
            return k
    return None


def tu_crate(mau_tep, mau_re):
    """Đọc một hằng số trong mã nguồn crate; None nếu máy không có crate."""
    goc = os.environ.get("CARGO_HOME", os.path.expanduser("~/.cargo"))
    tep = glob.glob(os.path.join(goc, "registry/src/*", mau_tep))
    if not tep:
        return None
    m = re.search(mau_re, open(tep[0], encoding="utf-8").read())
    return m.group(1) if m else None


def doi_chieu(ten, gia_tri, mau_tep, mau_re, co_so=10):
    s = tu_crate(mau_tep, mau_re)
    if s is None:
        print("  %s: không có crate để đối chiếu" % ten)
        return
    assert int(s, co_so) == gia_tri, ten
    print("  %s khớp crate %s" % (ten, mau_tep.split("/")[0]))


def in_duong_cong(ten, u, p, r, t, b, k_ky_vong):
    so = bac_nhom_j0(p, b, t)
    h, du = divmod(so, r)
    k = bac_nhung(p, r)
    assert du == 0 and miller_rabin(p) and miller_rabin(r) and k == k_ky_vong
    print("=== %s: y^2 = x^3 %s %d" % (ten, "-" if b < 0 else "+", abs(b)))
    if u is not None:
        print("u = %s%s (%d bit, trọng số Hamming %d)"
              % ("-" if u < 0 else "", hex(abs(u)), abs(u).bit_length(), bin(abs(u)).count("1")))
    print("p =", tu_64_bit(p), "(%d bit, 2-adicity %d)" % (p.bit_length(), two_adicity(p - 1)))
    print("r =", tu_64_bit(r), "(%d bit, 2-adicity %d)" % (r.bit_length(), two_adicity(r - 1)))
    print("h = %s (%d bit)" % (tu_64_bit(h), h.bit_length()))
    print("bậc nhúng k %s; p, r qua Miller-Rabin 40 cơ sở; #E = h*r đã kiểm trên 3 điểm"
          % ("= %d" % k if k else "> 10000"))
    return so


# --- BN254, BLS12-381: suy từ u ---
p, r, t = bn(U_BN254)
in_duong_cong("BN254", U_BN254, p, r, t, 3, 12)
doi_chieu("p", p, "ark-bn254-0.6.0/src/fields/fq.rs", r'modulus = "(\d+)"')
doi_chieu("r", r, "ark-bn254-0.6.0/src/fields/fr.rs", r'modulus = "(\d+)"')
doi_chieu("b", 3, "ark-bn254-0.6.0/src/curves/g1.rs", r'COEFF_B: Fq = MontFp!\("(\d+)"\)')

p, r, t = bls12(U_BLS12_381)
in_duong_cong("BLS12-381", U_BLS12_381, p, r, t, 4, 12)
doi_chieu("p", p, "ark-bls12-381-0.6.0/src/fields/fq.rs", r'modulus = "(\d+)"')
doi_chieu("r", r, "ark-bls12-381-0.6.0/src/fields/fr.rs", r'modulus = "(\d+)"')
doi_chieu("b", 4, "ark-bls12-381-0.6.0/src/curves/g1.rs", r'COEFF_B: Fq = MontFp!\("(\d+)"\)')

# --- BLS12-377: giá trị công bố phải khớp đa thức tại u ---
p, r, t = bls12(U_377)
assert (p, r) == (P_377, R_377) and U_377 % 3 == 1
in_duong_cong("BLS12-377", U_377, p, r, t, 1, 12)
print("  p, r công bố = p(u), r(u): True; u ≡ 1 (mod 3): True")

# --- BW6-761: dựng từ đa thức Bảng 4, trường vô hướng = trường cơ sở BLS12-377 ---
x = U_377
r6, du_r = divmod(x**6 - 2 * x**5 + 2 * x**3 + x + 1, 3)
t6 = x**5 - 3 * x**4 + 3 * x**3 - x + 3 + 13 * r6
y6, du_y = divmod(x**5 - 3 * x**4 + 3 * x**3 - x + 3, 3)
y6 += 9 * r6
q6, du_q = divmod(t6 * t6 + 3 * y6 * y6, 4)
q6_bang, du_b = divmod(103 * x**12 - 379 * x**11 + 250 * x**10 + 691 * x**9
                       - 911 * x**8 - 79 * x**7 + 623 * x**6 - 640 * x**5
                       + 274 * x**4 + 763 * x**3 + 73 * x**2 + 254 * x + 229, 9)
assert not (du_r or du_y or du_q or du_b)
assert q6 == q6_bang == Q_BW6 and r6 == P_377
in_duong_cong("BW6-761", None, q6, r6, t6, -1, 6)
so_xoan = bac_nhom_j0(q6, 4, t6)  # tập ứng viên không phụ thuộc dấu t
assert so_xoan % r6 == 0
print("  q = (t^2 + 3y^2)/4 = q(u) ở Bảng 4 = giá trị Bảng 3: True")
print("  r(BW6) = p(BLS12-377): True; twist y^2 = x^3 + 4 có r | #E': True")

# --- Pallas/Vesta: chu trình ---
assert P_PASTA == 2**254 + 45560315531419706090280762371685220353
assert Q_PASTA == 2**254 + 45560315531506369815346746415080538113
so_p = in_duong_cong("Pallas", None, P_PASTA, Q_PASTA, P_PASTA + 1 - Q_PASTA, 5, None)
so_q = in_duong_cong("Vesta", None, Q_PASTA, P_PASTA, Q_PASTA + 1 - P_PASTA, 5, None)
assert so_p == Q_PASTA and so_q == P_PASTA
print("  chu trình: #Pallas = q, #Vesta = p: True; q - p dài %d bit" % (Q_PASTA - P_PASTA).bit_length())
doi_chieu("p", P_PASTA, "halo2curves-0.9.0/src/pasta/fp.rs", r'modulus = "([0-9a-f]+)"', 16)
doi_chieu("q", Q_PASTA, "halo2curves-0.9.0/src/pasta/fq.rs", r'modulus = "([0-9a-f]+)"', 16)

# --- Số ứng viên chia hết cho r: một với năm đường cong, hai với BW6-761 (E và twist) ---
def so_ung_vien_chia_r(p, r, t):
    y = math.isqrt((4 * p - t * t) // 3)
    v1, v2 = (t + 3 * y) // 2, (t - 3 * y) // 2
    return sum((p + 1 - v) % r == 0 for v in {t, -t, v1, v2, -v1, -v2})


for bo in (bn(U_BN254), bls12(U_BLS12_381), bls12(U_377),
           (P_PASTA, Q_PASTA, P_PASTA + 1 - Q_PASTA), (Q_PASTA, P_PASTA, Q_PASTA + 1 - P_PASTA)):
    assert so_ung_vien_chia_r(*bo) == 1
assert so_ung_vien_chia_r(q6, r6, t6) == 2

# --- Đường cong toy: đếm trực tiếp, độc lập với snarklib ---
R, Q, A = 449, 5387, 1
dem = 1 + sum(1 + (0 if (v := (x**3 + x) % Q) == 0 else
                   (1 if pow(v, (Q - 1) // 2, Q) == 1 else -1)) for x in range(Q))
P_toy = (2, 1603)
G_toy = nhan(12, P_toy, A, Q)
assert dem == Q + 1 == 12 * R and Q % 4 == 3 and miller_rabin(Q) and miller_rabin(R)
assert G_toy == (2658, 1349) and nhan(R, G_toy, A, Q) is None and bac_nhung(Q, R) == 2
assert (R - 1) == 2**6 * 7 and all(pow(3, 448 // f, R) != 1 for f in (2, 7))
omega = {n: pow(3, 448 // n, R) for n in (8, 16, 32, 64)}
assert omega == {8: 122, 16: 349, 32: 221, 64: 391}
dong, loi = compute()
assert not loi and "e_G1_G1 = 1780 + 1636·i" in dong
print("=== Toy: y^2 = x^3 + x trên F_%d" % Q)
print("#E = %d = 12*%d (vết 0, siêu kỳ dị); bậc nhúng 2; g = 3" % (dem, R))
print("omega_8, 16, 32, 64 =", *omega.values())
print("P = %s, G1 = 12P = %s, e(G1, G1) = 1780 + 1636 i, (q^2-1)/r = %d"
      % (P_toy, G_toy, (Q * Q - 1) // R))
print("Tất cả phép kiểm đều đạt.")
