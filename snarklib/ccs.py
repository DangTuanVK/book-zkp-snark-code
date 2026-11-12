"""Customizable Constraint System (CCS) theo Setty–Thaler–Wahby.

Cấu trúc: t ma trận M_0..M_{t-1} cỡ m × n, q đa tập S_0..S_{q-1} (chỉ số ma trận),
q hằng số c_0..c_{q-1}. Vector z = (w, 1, x) theo bài báo gốc; để thống nhất
với r1cs.py, ở đây ta dùng z = (1, x, w) — chỉ là hoán vị cột, không đổi bản chất.
Thỏa khi:
    Σ_i c_i · ∘_{j ∈ S_i} (M_j z) = 0     (∘: tích Hadamard)
"""

from .r1cs import matvec


class CCS:
    def __init__(self, F, matrices, multisets, constants, n_public):
        self.F = F
        self.M = matrices
        self.S = [list(s) for s in multisets]
        self.c = [F(v) for v in constants]
        self.n_public = n_public
        self.m = len(matrices[0])
        self.n = len(matrices[0][0])
        assert len(self.S) == len(self.c)
        for M in matrices:
            assert len(M) == self.m and all(len(row) == self.n for row in M)

    def z(self, x, w):
        if len(x) != self.n_public or 1 + len(x) + len(w) != self.n:
            raise ValueError("độ dài x hoặc w không khớp")
        return [self.F(1)] + [self.F(v) for v in x] + [self.F(v) for v in w]

    def residuals(self, x, w):
        F = self.F
        z = self.z(x, w)
        Mz = [matvec(F, M, z) for M in self.M]
        out = [F(0)] * self.m
        for ci, Si in zip(self.c, self.S):
            prod = [F(1)] * self.m
            for j in Si:
                prod = [p * v for p, v in zip(prod, Mz[j])]
            out = [o + ci * p for o, p in zip(out, prod)]
        return out

    def is_satisfied(self, x, w):
        return all(v == 0 for v in self.residuals(x, w))

    def degree(self):
        """Bậc của CCS: max |S_i|."""
        return max(len(s) for s in self.S)


def from_r1cs(r):
    """R1CS -> CCS: M = (A, B, C), S = ({0, 1}, {2}), c = (1, -1)."""
    return CCS(r.F, [r.A, r.B, r.C], [[0, 1], [2]], [1, -1], r.n_public)
