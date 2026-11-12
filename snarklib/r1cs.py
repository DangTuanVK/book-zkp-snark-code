"""Hệ ràng buộc bậc hai R1CS: (A z) ∘ (B z) = C z, với z = (1, x, w).

- A, B, C: ma trận m × n (danh sách m hàng, mỗi hàng n số), n = 1 + len(x) + len(w).
- x: đầu vào công khai (instance), w: witness.
- ∘: tích Hadamard (nhân từng thành phần).
"""


def matvec(F, M, z):
    """Tích ma trận–vector trên F."""
    return [sum((F(a) * zi for a, zi in zip(row, z)), F(0)) for row in M]


class R1CS:
    """Cấu trúc R1CS tối giản."""

    def __init__(self, F, A, B, C, n_public):
        self.F = F
        self.A, self.B, self.C = A, B, C
        self.m = len(A)
        self.n = len(A[0]) if A else 0
        self.n_public = n_public          # số phần tử của x
        for M in (A, B, C):
            assert len(M) == self.m and all(len(row) == self.n for row in M)

    def z(self, x, w):
        """Ghép z = (1, x, w)."""
        if len(x) != self.n_public or 1 + len(x) + len(w) != self.n:
            raise ValueError("độ dài x hoặc w không khớp")
        return [self.F(1)] + [self.F(v) for v in x] + [self.F(v) for v in w]

    def residuals(self, x, w):
        """Vector (Az)∘(Bz) - Cz; toàn 0 khi thỏa."""
        z = self.z(x, w)
        Az = matvec(self.F, self.A, z)
        Bz = matvec(self.F, self.B, z)
        Cz = matvec(self.F, self.C, z)
        return [a * b - c for a, b, c in zip(Az, Bz, Cz)]

    def is_satisfied(self, x, w):
        return all(v == 0 for v in self.residuals(x, w))


def factor_r1cs(F):
    """Ví dụ một ràng buộc a·b = N: z = (1, N, a, b), chưa có kiểm khoảng."""
    A = [[0, 0, 1, 0]]
    B = [[0, 0, 0, 1]]
    C = [[0, 1, 0, 0]]
    return R1CS(F, A, B, C, n_public=1)
