"""KZG đa tuyến tính của Papamanthou–Shi–Tamassia (PST13), theo trình bày ở
[Thaler, Proofs, Arguments, and Zero-Knowledge, §15.3].

Đa thức đa tuyến tính l biến q cho bởi bảng 2^l giá trị trên hypercube (quy ước
big-endian của snarklib.mle: x_1 là bit cao). Cửa sập là véc-tơ τ = (τ_1..τ_l).

    SRS:     tables[j] = [ eq(b, (τ_{j+1},...,τ_l))·G1  với mọi b ∈ {0,1}^(l-j) ],
             j = 0..l  (tables[0] dùng để commit, tables[j] để cam kết thương w_j)
             g2 = [ [1]_2, [τ_1]_2, ..., [τ_l]_2 ]
    Commit:  C = Σ_b q(b)·[eq(b, τ)]_1 = [q(τ)]_1
    Open(z): q(X) - v = Σ_i (X_i - z_i)·w_i(X_{i+1},...,X_l)       (Fact 15.1)
             w_i = MLE của bảng hiệu T_hi - T_lo sau khi đã fold z_1..z_{i-1}
             π = (W_1..W_l), W_i = [w_i(τ_{i+1..l})]_1
    Verify:  e(C - [v]_1, [1]_2) = Π_i e(W_i, [τ_i]_2 - z_i·[1]_2)

Toy: pairing đối xứng (G2 = G1); tham số G1, G2 tách bạch để viết dạng Type 3.
"""

from .field import GF
from .curve import toy_generator, R_TOY
from .mle import eq_table, evaluate, fold, num_vars
from . import opcount as oc

F = GF(R_TOY)


class PSTSRS:
    def __init__(self, tables, g2):
        self.tables, self.g2 = tables, g2

    @property
    def num_vars(self):
        return len(self.g2) - 1

    @property
    def G1(self):
        return self.tables[-1][0]


def setup(l, taus, G1=None, G2=None):
    """SRS cho đa thức l biến từ cửa sập taus = (τ_1..τ_l) (không đếm opcount)."""
    if len(taus) != l:
        raise ValueError("cần đúng l giá trị τ")
    if G1 is None:
        G1 = toy_generator()[2]
    if G2 is None:
        G2 = G1                       # toy: Type 1
    taus = [F(t) for t in taus]
    tables = []
    for j in range(l + 1):
        tables.append([int(e) * G1 for e in eq_table(F, taus[j:])])
    g2 = [G2] + [int(t) * G2 for t in taus]
    return PSTSRS(tables, g2)


def commit(srs, table):
    """C = [q(τ)]_1: một MSM cỡ 2^l trên cơ sở Lagrange đã mã hóa."""
    if len(table) != len(srs.tables[0]):
        raise ValueError("bảng cần 2^l phần tử")
    return oc.msm(srs.tables[0], [F(v) for v in table])


def open(srs, table, z):
    """Trả về (v = q(z), [W_1..W_l]); W_i cam kết thương w_i (MSM cỡ 2^(l-i))."""
    l = num_vars(table)
    if len(z) != l:
        raise ValueError("điểm mở phải có l tọa độ")
    T = [F(v) for v in table]
    Ws = []
    for i in range(l):
        half = len(T) // 2
        diff = [T[j + half] - T[j] for j in range(half)]   # bảng của w_{i+1}
        Ws.append(oc.msm(srs.tables[i + 1], diff))
        T = fold(F, T, z[i])                               # cố định x_{i+1} = z_{i+1}
    v = T[0]
    assert v == evaluate(F, table, z)
    return v, Ws


def quotients(table, z):
    """Các bảng thương w_1..w_l (để in ví dụ số và kiểm Fact 15.1)."""
    T = [F(v) for v in table]
    out = []
    for i in range(num_vars(table)):
        half = len(T) // 2
        out.append([T[j + half] - T[j] for j in range(half)])
        T = fold(F, T, z[i])
    return out


def verify(srs, C, z, v, Ws):
    """Kiểm e(C - [v]_1, [1]_2) · Π_i e(-W_i, [τ_i]_2 - z_i[1]_2) = 1 (l + 1 pairing)."""
    l = srs.num_vars
    if len(z) != l or len(Ws) != l:
        return False
    pairs = [(oc.msm([C, srs.G1], [1, -F(v)]), srs.g2[0])]
    for i in range(l):
        Qi = oc.msm([srs.g2[i + 1], srs.g2[0]], [1, -F(z[i])], group="G2")
        pairs.append((-Ws[i], Qi))
    return oc.pairing_check(pairs)
