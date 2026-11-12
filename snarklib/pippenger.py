"""MSM (multi-scalar multiplication) Σ s_i·P_i: Pippenger (phương pháp xô) và ngây thơ.

Cả hai thuật toán chạy trên một "nhóm" trừu tượng cung cấp add(a, b) và dbl(a);
phần tử trung hòa *cấu trúc* (xô rỗng, bộ tích lũy chưa khởi tạo) biểu diễn
bằng None và không tính phí. Nhờ vậy số phép toán đếm được chỉ phụ thuộc chữ
số của vô hướng, không phụ thuộc đường cong:
    NhomDiem(E)   tính thật trên đường cong toy (curve.Point) -> kiểm kết quả
    NhomDem()     không tính gì, chỉ đếm -> quét n, c lớn nhanh

Pippenger với cửa sổ c bit, vô hướng b bit, W = ceil(b/c) cửa sổ:
    với mỗi cửa sổ k (từ cao xuống thấp):
        acc = c lần nhân đôi acc                      (bỏ qua ở cửa sổ đầu)
        xô B_1..B_{2^c - 1}: B_d += P_i với mọi i có chữ số d_{i,k} = d != 0
        tổng có trọng số Σ d·B_d bằng tổng chạy:
            S = T = O; với d = 2^c - 1 xuống 1: S += B_d; T += S
        acc += T
Chi phí lý thuyết (cận trên, mỗi phép trên xô/tổng chạy đều tính):
    cộng   ≈ W·(n + 2^{c+1})          nhân đôi ≈ c·(W - 1) ≈ b
Ngây thơ: mỗi điểm double-and-add (b-1 nhân đôi, wt(s_i)-1 phép cộng)
cộng thêm n - 1 phép cộng gộp; trung bình ≈ n·(b - 1) nhân đôi + n·b/2 phép cộng.

Đường cong toy là pairing đối xứng (Type 1, G1 = G2); MSM ở đây chỉ dùng G1.
Thuật toán: Pippenger 1976 (SIAM J. Comput.); dạng "xô + tổng chạy" theo các
cài đặt hiện đại (vd. arkworks ark-ec VariableBaseMSM).
"""


class DemPhepToan:
    """Bộ đếm phép cộng điểm (add) và nhân đôi (dbl)."""

    def __init__(self):
        self.add = 0
        self.dbl = 0

    def tong(self):
        return self.add + self.dbl

    def __repr__(self):
        return "cộng = %d, nhân đôi = %d" % (self.add, self.dbl)


class NhomDiem:
    """Phép toán thật trên điểm của curve.Curve (affine), có đếm."""

    def __init__(self):
        self.dem = DemPhepToan()

    def add(self, a, b):
        if a is None:
            return b
        if b is None:
            return a
        self.dem.add += 1
        return a + b

    def dbl(self, a):
        if a is None:
            return None
        self.dem.dbl += 1
        return a + a


class NhomDem:
    """Chỉ đếm: phần tử là token True, None là trung hòa cấu trúc."""

    def __init__(self):
        self.dem = DemPhepToan()

    def add(self, a, b):
        if a is None:
            return b
        if b is None:
            return a
        self.dem.add += 1
        return True

    def dbl(self, a):
        if a is None:
            return None
        self.dem.dbl += 1
        return True


def _ket_qua(acc, points):
    """Đổi trung hòa cấu trúc None thành điểm vô cực O khi tính trên điểm thật."""
    if acc is None and points and hasattr(points[0], "E"):
        return points[0].E.infinity()
    return acc


def msm_ngay_tho(G, points, scalars):
    """Σ s_i·P_i: double-and-add từng điểm (bit cao xuống thấp) rồi cộng dồn."""
    acc = None
    for P, s in zip(points, scalars):
        R = None
        for bit in bin(s)[2:] if s > 0 else "":
            R = G.dbl(R)
            if bit == "1":
                R = G.add(R, P)
        acc = G.add(acc, R)
    return _ket_qua(acc, points)


def msm_pippenger(G, points, scalars, c, bits=None):
    """Σ s_i·P_i bằng phương pháp xô với cửa sổ c bit (vô hướng không âm, <= bits bit)."""
    if bits is None:
        bits = max(s.bit_length() for s in scalars) or 1
    W = -(-bits // c)
    mask = (1 << c) - 1
    acc = None
    for k in range(W - 1, -1, -1):
        for _ in range(c):
            acc = G.dbl(acc)
        buckets = [None] * (1 << c)          # buckets[0] không dùng
        for P, s in zip(points, scalars):
            d = (s >> (c * k)) & mask
            if d:
                buckets[d] = G.add(buckets[d], P)
        S = T = None
        for d in range(mask, 0, -1):          # tổng chạy: T = Σ d·B_d
            S = G.add(S, buckets[d])
            T = G.add(T, S)
        acc = G.add(acc, T)
    return _ket_qua(acc, points)


def chi_phi_ly_thuyet(n, b, c):
    """(số phép cộng, số nhân đôi) ước lượng: W·(n + 2^{c+1}), c·(W - 1)."""
    W = -(-b // c)
    return W * (n + (1 << (c + 1))), c * (W - 1)


def c_toi_uu_ly_thuyet(n, b, c_max=24):
    """c cực tiểu hóa W·(n + 2^{c+1}) + c·(W - 1)."""
    return min(range(1, c_max + 1), key=lambda c: sum(chi_phi_ly_thuyet(n, b, c)))


def dem_pippenger(scalars, c, bits):
    """Đếm phép toán của Pippenger (không tính điểm)."""
    G = NhomDem()
    msm_pippenger(G, [True] * len(scalars), scalars, c, bits)
    return G.dem


def dem_ngay_tho(scalars):
    G = NhomDem()
    msm_ngay_tho(G, [True] * len(scalars), scalars)
    return G.dem
