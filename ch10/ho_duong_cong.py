"""Hàm phụ trợ Ch.10 (không in nguyên văn trong sách): kiểm nguyên tố
Miller–Rabin, đa thức tham số BN/BLS12, bậc nhúng, 2-adicity.

Giá trị u đối chiếu với mã nguồn arkworks 0.6.0 (registry cargo):
  BN254:     ark-bn254-0.6.0/src/curves/mod.rs:18-20
             X = [4965661367192848881], X_IS_NEGATIVE = false
  BLS12-381: ark-bls12-381-0.6.0/src/curves/mod.rs:25-26
             X = [0xd201000000010000], X_IS_NEGATIVE = true
p, r đối chiếu: ark-bn254-0.6.0/src/fields/fq.rs:4, fr.rs:4;
                ark-bls12-381-0.6.0/src/fields/fq.rs:4, fr.rs:4.
"""

# ---- Hằng số chép từ arkworks 0.6.0 (xem docstring) ----
U_BN254 = 4965661367192848881
U_BLS12_381 = -0xd201000000010000
ARK_P_BN254 = 21888242871839275222246405745257275088696311157297823662689037894645226208583
ARK_R_BN254 = 21888242871839275222246405745257275088548364400416034343698204186575808495617
ARK_P_BLS12_381 = int("40024095552216673934177898257359041565568828199390078853320581361"
                      "24031650490837864442687629129015664037894272559787")
ARK_R_BLS12_381 = 52435875175126190479447740508185965837690552500527637822603658699938581184513
# ark-bls12-381-0.6.0/src/curves/g1.rs:32-33: COFACTOR = (x-1)^2/3
ARK_H1_BLS12_381 = 0x396c8c005555e1568c00aaab0000aaab
# ark-bn254-0.6.0/src/curves/mod.rs:21-25: ATE_LOOP_COUNT (NAF của 6u+2, bit thấp trước)
ARK_ATE_LOOP_BN254 = [
    0, 0, 0, 1, 0, 1, 0, -1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, -1, 0, 0, 0, 1, 0, -1, 0, 0, 0,
    0, -1, 0, 0, 1, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0, 1, 0, -1, 0, 0, 0, -1, 0,
    -1, 0, 0, 0, 1, 0, 1, 1,
]

# 40 số nguyên tố đầu tiên làm cơ sở Miller–Rabin
CO_SO_MR = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
            73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151,
            157, 163, 167, 173]


def miller_rabin(n, bases=CO_SO_MR):
    """Miller–Rabin với tập cơ sở cố định (tất định, không dùng số ngẫu nhiên).

    Với n < 3.3·10^24 thì 13 cơ sở đầu đã đủ chứng minh nguyên tố; với số 254-381 bit,
    40 cơ sở cho xác suất hợp số lọt qua <= 4^-40 nếu n không được dựng đối kháng.
    """
    if n < 2:
        return False
    for q in bases:
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in bases:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def bn(u):
    """Họ BN (Barreto–Naehrig 2005): p(u), r(u), vết t(u) = 6u^2 + 1."""
    p = 36 * u**4 + 36 * u**3 + 24 * u**2 + 6 * u + 1
    r = 36 * u**4 + 36 * u**3 + 18 * u**2 + 6 * u + 1
    return p, r, 6 * u**2 + 1


def bls12(u):
    """Họ BLS12 (Barreto–Lynn–Scott 2002): p(u) = (u-1)^2 (u^4-u^2+1)/3 + u, r(u) = u^4-u^2+1, t = u+1."""
    num = (u - 1) ** 2 * (u**4 - u**2 + 1)
    assert num % 3 == 0
    return num // 3 + u, u**4 - u**2 + 1, u + 1


def bac_nhung(p, r, k_max=50):
    """Bậc nhúng = bậc của p trong F_r^*: k nhỏ nhất với r | p^k - 1."""
    x = 1
    for k in range(1, k_max + 1):
        x = x * p % r
        if x == 1:
            return k
    return None


def two_adicity(n):
    s = 0
    while n % 2 == 0:
        n //= 2
        s += 1
    return s


def naf_gia_tri(digits):
    """Giá trị của dãy chữ số NAF (bit thấp trước)."""
    return sum(d << i for i, d in enumerate(digits))
