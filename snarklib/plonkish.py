"""Số học hóa Plonkish (kiểu PLONK gốc, 3 cột a, b, c).

Mỗi hàng j thỏa cổng:
    q_L[j]·a[j] + q_R[j]·b[j] + q_O[j]·c[j] + q_M[j]·a[j]·b[j] + q_C[j] + PI[j] = 0
Ràng buộc sao chép: danh sách các cặp ô ((cột, hàng), (cột, hàng)) phải bằng nhau,
cột là một trong "a", "b", "c". Theo quy ước PLONK, đầu vào công khai x_i đặt
ở hàng i với q_L = 1, a = x_i và PI[i] = -x_i.
"""

SELECTORS = ("qL", "qR", "qO", "qM", "qC")
COLUMNS = ("a", "b", "c")


class PlonkishCircuit:
    """Mạch Plonkish: selector cố định + ràng buộc sao chép."""

    def __init__(self, F, rows, copies=(), n_public=0):
        """rows: danh sách dict {qL, qR, qO, qM, qC} (thiếu khóa = 0)."""
        self.F = F
        self.n = len(rows)
        self.sel = {s: [F(r.get(s, 0)) for r in rows] for s in SELECTORS}
        self.copies = list(copies)
        self.n_public = n_public

    def public_input_column(self, x):
        """Cột PI: PI[i] = -x_i cho i < len(x), còn lại 0."""
        if len(x) != self.n_public:
            raise ValueError("số đầu vào công khai không khớp")
        return [-self.F(x[i]) if i < len(x) else self.F(0) for i in range(self.n)]

    def gate_residuals(self, a, b, c, x=()):
        """Giá trị vế trái của cổng ở từng hàng."""
        F, s = self.F, self.sel
        PI = self.public_input_column(list(x))
        a, b, c = [F(v) for v in a], [F(v) for v in b], [F(v) for v in c]
        return [s["qL"][j] * a[j] + s["qR"][j] * b[j] + s["qO"][j] * c[j]
                + s["qM"][j] * a[j] * b[j] + s["qC"][j] + PI[j]
                for j in range(self.n)]

    def copy_violations(self, a, b, c):
        """Danh sách các cặp sao chép bị vi phạm."""
        col = {"a": [self.F(v) for v in a], "b": [self.F(v) for v in b],
               "c": [self.F(v) for v in c]}
        bad = []
        for (c1, r1), (c2, r2) in self.copies:
            if col[c1][r1] != col[c2][r2]:
                bad.append(((c1, r1), (c2, r2)))
        return bad

    def is_satisfied(self, a, b, c, x=()):
        if not (len(a) == len(b) == len(c) == self.n):
            raise ValueError("độ dài cột khác số hàng")
        gates_ok = all(v == 0 for v in self.gate_residuals(a, b, c, x))
        return gates_ok and not self.copy_violations(a, b, c)


def factor_plonkish(F):
    """Ví dụ: công khai N, chứng minh biết a, b với a·b = N.

    Hàng 0: cổng công khai  a0 - N = 0      (qL = 1, PI = -N)
    Hàng 1: cổng nhân       a1·b1 - c1 = 0  (qM = 1, qO = -1)
    Sao chép: c1 = a0.
    """
    rows = [{"qL": 1}, {"qM": 1, "qO": -1}]
    copies = [(("c", 1), ("a", 0))]
    return PlonkishCircuit(F, rows, copies, n_public=1)
