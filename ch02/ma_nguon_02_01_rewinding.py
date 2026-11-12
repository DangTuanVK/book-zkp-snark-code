"""Extractor tua lại cho giao thức Schnorr trên nhóm toy (Chương 2).
Đo xác suất trích xuất thành công theo xác suất thắng eps của đối thủ,
so với cận eps^2 - eps/N của bổ đề tua lại (N = 449 thách thức)."""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.field import GF
from snarklib.curve import toy_generator

r = 449
F = GF(r)
E, _, G = toy_generator()
w = 17                                    # witness: U = w·G1
U = G * w
bG = [G * i for i in range(r)]            # bảng bội của G1 và U
bU = [U * i for i in range(r)]


def chap_nhan(T, c, z):                   # Verifier: z·G1 = T + c·U ?
    return bG[z] == T + bU[c]


class DoiThu:
    """Prover gian lận: biết w nhưng chỉ trả lời đúng với tập thách thức S."""
    def __init__(self, kieu, m, rng):
        self.kieu, self.m, self.rng = kieu, m, rng

    def cam_ket(self):
        self.k = self.rng.randrange(r)
        if self.kieu == "deu":            # mọi trạng thái: đúng m thách thức
            self.S = set(random.Random(self.k).sample(range(r), self.m))
        else:                             # tất cả hoặc không: P[S = F] = m/r
            tot = self.rng.randrange(r) < self.m
            self.S = set(range(r)) if tot else set()
        return bG[self.k]

    def tra_loi(self, c):                 # trạng thái (k, S) giữ nguyên khi tua lại
        z = (self.k + c * w) % r
        return z if c in self.S else (z + 1) % r


def trich_xuat(A, rng):
    """Một lần chạy, tua lại về sau cam kết, chạy lần hai với thách thức mới."""
    T = A.cam_ket()
    c1 = rng.randrange(r); z1 = A.tra_loi(c1)
    c2 = rng.randrange(r); z2 = A.tra_loi(c2)
    if c1 != c2 and chap_nhan(T, c1, z1) and chap_nhan(T, c2, z2):
        w_hat = int(F(z1 - z2) / F(c1 - c2))
        assert bG[w_hat] == U             # witness trích được là đúng
        return w_hat
    return None


rng = random.Random(2026)
LAN = 4000
print("  m    eps   | deu: thang  trich   can  | tat-ca: thang  trich")
for m in (22, 45, 90, 135, 180, 225, 270, 315, 360, 405, 449):
    eps = m / r
    hang = []
    for kieu in ("deu", "tat-ca"):
        A = DoiThu(kieu, m, rng)
        thang = sum(chap_nhan(A.cam_ket(), c, A.tra_loi(c))
                    for c in (rng.randrange(r) for _ in range(LAN)))
        trich = sum(trich_xuat(A, rng) is not None for _ in range(LAN))
        hang.append((thang / LAN, trich / LAN))
    can = eps * eps - eps / r             # bổ đề tua lại
    assert hang[0][1] >= can - 0.03 and hang[1][1] >= can - 0.03
    print("%4d  %.3f  |   %.3f   %.3f  %.3f |    %.3f   %.3f"
          % (m, eps, hang[0][0], hang[0][1], can, hang[1][0], hang[1][1]))
