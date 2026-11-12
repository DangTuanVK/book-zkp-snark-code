"""Đường cong toy E: y^2 = x^3 + x trên F_q, q = 5387 (siêu kỳ dị).

#E(F_q) = q + 1 = 5388 = 12 · 449, nhóm con bậc r = 449 dùng cho pairing.
Hàm cũng nhận tham số tổng quát E: y^2 = x^3 + A x + B trên F_q.
"""

from .field import GF

Q_TOY = 5387
R_TOY = 449
COFACTOR_TOY = 12


class Curve:
    """Đường cong Weierstrass ngắn y^2 = x^3 + A x + B trên F_q."""

    def __init__(self, q=Q_TOY, A=1, B=0):
        self.F = GF(q)
        self.q = q
        self.A = self.F(A)
        self.B = self.F(B)
        if 4 * self.A ** 3 + 27 * self.B ** 2 == 0:
            raise ValueError("đường cong suy biến")

    def __repr__(self):
        return "E: y^2 = x^3 + %d·x + %d trên F_%d" % (int(self.A), int(self.B), self.q)

    def rhs(self, x):
        x = self.F(x)
        return x ** 3 + self.A * x + self.B

    def infinity(self):
        return Point(self, None, None)

    def point(self, x, y):
        P = Point(self, self.F(x), self.F(y))
        if not P.on_curve():
            raise ValueError("(%d, %d) không nằm trên đường cong" % (int(x), int(y)))
        return P

    def count_points(self):
        """Đếm #E(F_q) = 1 + Σ_x (1 + Legendre(x^3 + Ax + B)) (duyệt mọi x)."""
        total = 1  # điểm vô cực
        for x in range(self.q):
            total += 1 + self.rhs(x).legendre()
        return total

    def points(self):
        """Liệt kê mọi điểm affine (không gồm O) — chỉ dùng với q nhỏ."""
        for x in range(self.q):
            v = self.rhs(x)
            if v == 0:
                yield Point(self, self.F(x), self.F(0))
            elif v.legendre() == 1:
                y = v.sqrt()
                yield Point(self, self.F(x), y)
                yield Point(self, self.F(x), -y)

    def deterministic_generator(self, r=R_TOY, cofactor=COFACTOR_TOY):
        """Sinh G1 tất định cho nhóm con bậc r.

        Duyệt x = 0, 1, 2, ...: lấy x nhỏ nhất mà x^3+Ax+B là thặng dư bậc hai
        khác 0, y = căn nhỏ hơn, P = (x, y), G1 = cofactor·P; nếu G1 = O thì
        sang x kế tiếp. Trả về (P, G1).
        """
        for x in range(self.q):
            v = self.rhs(x)
            if v == 0 or v.legendre() != 1:
                continue
            P = Point(self, self.F(x), v.sqrt())
            G = cofactor * P
            if G.is_infinity():
                continue
            assert (r * G).is_infinity()
            return P, G
        raise ValueError("không tìm được điểm sinh")


class Point:
    """Điểm affine (x, y) trên Curve; x = y = None biểu diễn điểm vô cực O."""

    __slots__ = ("E", "x", "y")

    def __init__(self, E, x, y):
        self.E, self.x, self.y = E, x, y

    def is_infinity(self):
        return self.x is None

    def on_curve(self):
        if self.is_infinity():
            return True
        return self.y ** 2 == self.E.rhs(self.x)

    def __repr__(self):
        if self.is_infinity():
            return "O"
        return "(%d, %d)" % (int(self.x), int(self.y))

    def __eq__(self, o):
        if not isinstance(o, Point):
            return False
        if self.is_infinity() or o.is_infinity():
            return self.is_infinity() and o.is_infinity()
        return self.x == o.x and self.y == o.y

    def __ne__(self, o):
        return not self.__eq__(o)

    def __hash__(self):
        return hash(None) if self.is_infinity() else hash((int(self.x), int(self.y)))

    def __neg__(self):
        if self.is_infinity():
            return self
        return Point(self.E, self.x, -self.y)

    def __add__(self, o):
        """Cộng điểm affine (luật hợp dây cung–tiếp tuyến)."""
        if self.is_infinity():
            return o
        if o.is_infinity():
            return self
        if self.x == o.x and self.y == -o.y:
            return self.E.infinity()                      # P + (-P) = O
        if self == o:
            lam = (3 * self.x ** 2 + self.E.A) / (2 * self.y)  # tiếp tuyến
        else:
            lam = (o.y - self.y) / (o.x - self.x)              # dây cung
        x3 = lam ** 2 - self.x - o.x
        y3 = lam * (self.x - x3) - self.y
        return Point(self.E, x3, y3)

    def __sub__(self, o):
        return self + (-o)

    def __mul__(self, k):
        """Nhân vô hướng double-and-add (duyệt bit từ cao xuống thấp)."""
        k = int(k)
        if k < 0:
            return (-self) * (-k)
        R = self.E.infinity()
        for bit in bin(k)[2:]:
            R = R + R
            if bit == "1":
                R = R + self
        return R

    __rmul__ = __mul__

    def order(self, candidates=None):
        """Bậc của điểm: ước nhỏ nhất d của #E (hoặc của candidates) với d·P = O."""
        n = candidates if candidates is not None else self.E.count_points()
        for d in range(1, n + 1):
            if n % d == 0 and (d * self).is_infinity():
                return d
        raise ValueError("không tìm được bậc")


def toy_curve():
    """Đường cong toy của sách: y^2 = x^3 + x trên F_5387."""
    return Curve(Q_TOY, 1, 0)


def toy_generator():
    """Trả về (E, P, G1) theo quy tắc sinh tất định."""
    E = toy_curve()
    P, G = E.deterministic_generator()
    return E, P, G
