"""Cam kết Pedersen véc-tơ trên nhóm con bậc r = 449 của đường cong toy.

    C = Com(a; ρ) = Σ_i a_i·G_i + ρ·H        (ký hiệu cộng)

Véc-tơ cơ sở G_0, ..., G_{n-1}, H, U được sinh TẤT ĐỊNH bằng băm-vào-đường-cong
kiểu thử-và-tăng (try-and-increment), để "không ai biết quan hệ logarit" giữa
chúng (không có cửa sau, không cần trusted setup):

    với ctr = 0, 1, 2, ...:
        x = SHA256(DOMAIN ‖ nhãn ‖ chỉ số ‖ ctr)  mod q
        nếu x^3 + x là thặng dư bậc hai khác 0:
            y = căn NHỎ HƠN (quy ước sqrt của snarklib), P = (x, y)
            G = 12·P (nhân cofactor để rơi vào nhóm con bậc 449)
            nếu G ≠ O: trả về G

CẢNH BÁO: nhóm chỉ có 449 phần tử nên logarit rời rạc giải được tức thì bằng
vét cạn (xem dlog_bruteforce) — tính binding ở đây chỉ để MINH HỌA cơ chế.
Hệ thật dùng nhóm bậc ~2^255 (Pallas/Vesta, Ristretto, BLS12-381 G1).
Lưu ý: Pedersen chỉ cần một nhóm (không cần pairing); đường cong toy là
pairing Type 1 nhưng ở đây ta chỉ dùng phép cộng điểm.
"""

import hashlib

from .curve import toy_curve, R_TOY, COFACTOR_TOY

DOMAIN = b"snarklib/pedersen/v1"
E = toy_curve()


class OpCount:
    """Bộ đếm phép toán nhóm cục bộ: smul = nhân vô hướng, add = cộng điểm."""

    def __init__(self):
        self.smul = 0
        self.add = 0

    def as_dict(self):
        return {"smul": self.smul, "add": self.add}

    def __repr__(self):
        return "OpCount(smul=%d, add=%d)" % (self.smul, self.add)


def smul(P, k, cnt=None):
    """k·P (k đưa về [0, r)); đếm vào cnt nếu có."""
    if cnt is not None:
        cnt.smul += 1
    return int(k) % R_TOY * P


def padd(P, Q, cnt=None):
    if cnt is not None:
        cnt.add += 1
    return P + Q


def msm(scalars, points, cnt=None):
    """Σ k_i·P_i bằng cách ngây thơ (n nhân vô hướng + n-1 cộng điểm)."""
    if len(scalars) != len(points):
        raise ValueError("độ dài khác nhau")
    acc = None
    for k, P in zip(scalars, points):
        t = smul(P, k, cnt)
        acc = t if acc is None else padd(acc, t, cnt)
    return acc if acc is not None else E.infinity()


def hash_to_curve(label, index):
    """Băm (nhãn, chỉ số) vào nhóm con bậc 449. Trả về (G, ctr đã dùng)."""
    if isinstance(label, str):
        label = label.encode("utf-8")
    for ctr in range(1 << 16):
        h = hashlib.sha256(DOMAIN + len(label).to_bytes(4, "big") + label
                           + index.to_bytes(4, "big") + ctr.to_bytes(4, "big")).digest()
        x = int.from_bytes(h, "big") % E.q
        v = E.rhs(x)
        if v == 0 or v.legendre() != 1:
            continue
        P = E.point(x, v.sqrt())
        G = COFACTOR_TOY * P
        if G.is_infinity():
            continue
        assert (R_TOY * G).is_infinity()
        return G, ctr
    raise ValueError("không tìm được điểm")


def generators(label, n):
    """n điểm cơ sở độc lập G_0..G_{n-1} sinh từ nhãn."""
    return [hash_to_curve(label, i)[0] for i in range(n)]


class PedersenParams:
    """Tham số công khai: G (véc-tơ n điểm), H (cho ngẫu nhiên hóa), U (cho tích vô hướng)."""

    def __init__(self, n, label=b"snarklib/pedersen"):
        self.n = n
        self.label = label
        self.G = generators(label + b"/G", n)
        self.H = hash_to_curve(label + b"/H", 0)[0]
        self.U = hash_to_curve(label + b"/U", 0)[0]


def commit(G, a, H=None, rho=0, cnt=None):
    """Com(a; ρ) = <a, G> + ρ·H. Với H = None: cam kết không ẩn (non-hiding)."""
    if len(a) > len(G):
        raise ValueError("véc-tơ dài hơn số điểm cơ sở")
    C = msm(a, G[:len(a)], cnt)
    if H is not None and int(rho) % R_TOY:
        C = padd(C, smul(H, rho, cnt), cnt)
    return C


def dlog_bruteforce(base, P):
    """Tìm k ∈ [0, r) với k·base = P bằng vét cạn — CHỈ khả thi vì r = 449."""
    acc = E.infinity()
    for k in range(R_TOY):
        if acc == P:
            return k
        acc = acc + base
    raise ValueError("P không thuộc nhóm sinh bởi base")
