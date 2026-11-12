"""Đa thức một biến trên F_p.

Hệ số lưu từ bậc thấp đến bậc cao: Poly(F, [c0, c1, c2]) = c0 + c1·X + c2·X^2.
"""


class Poly:
    """Đa thức một biến với hệ số trong trường F (lớp do GF(p) tạo ra)."""

    def __init__(self, F, coeffs):
        self.F = F
        c = [F(x) for x in coeffs]
        while c and c[-1] == 0:      # bỏ hệ số 0 ở bậc cao
            c.pop()
        self.coeffs = c

    # --- tiện ích ---
    @classmethod
    def zero(cls, F):
        return cls(F, [])

    @classmethod
    def const(cls, F, c):
        return cls(F, [c])

    @classmethod
    def X(cls, F):
        return cls(F, [0, 1])

    def degree(self):
        """Bậc của đa thức; đa thức 0 có bậc -1 theo quy ước."""
        return len(self.coeffs) - 1

    def is_zero(self):
        return not self.coeffs

    def __repr__(self):
        if not self.coeffs:
            return "0"
        terms = []
        for k, c in enumerate(self.coeffs):
            if c == 0:
                continue
            if k == 0:
                terms.append("%d" % int(c))
            elif k == 1:
                terms.append("%d·X" % int(c))
            else:
                terms.append("%d·X^%d" % (int(c), k))
        return " + ".join(terms)

    def __eq__(self, o):
        if not isinstance(o, Poly):
            o = Poly(self.F, [o])
        return [int(x) for x in self.coeffs] == [int(x) for x in o.coeffs]

    def _lift(self, o):
        return o if isinstance(o, Poly) else Poly(self.F, [o])

    # --- số học ---
    def __add__(self, o):
        o = self._lift(o)
        n = max(len(self.coeffs), len(o.coeffs))
        a = self.coeffs + [self.F(0)] * (n - len(self.coeffs))
        b = o.coeffs + [self.F(0)] * (n - len(o.coeffs))
        return Poly(self.F, [x + y for x, y in zip(a, b)])

    __radd__ = __add__

    def __neg__(self):
        return Poly(self.F, [-c for c in self.coeffs])

    def __sub__(self, o):
        return self + (-self._lift(o))

    def __rsub__(self, o):
        return self._lift(o) - self

    def __mul__(self, o):
        o = self._lift(o)
        if self.is_zero() or o.is_zero():
            return Poly.zero(self.F)
        out = [self.F(0)] * (len(self.coeffs) + len(o.coeffs) - 1)
        for i, a in enumerate(self.coeffs):
            for j, b in enumerate(o.coeffs):
                out[i + j] = out[i + j] + a * b
        return Poly(self.F, out)

    __rmul__ = __mul__

    def __pow__(self, e):
        result, base = Poly.const(self.F, 1), self
        while e:
            if e & 1:
                result = result * base
            base = base * base
            e >>= 1
        return result

    def divmod(self, d):
        """Chia có dư: trả về (thương, dư) với self = thương·d + dư, deg(dư) < deg(d)."""
        d = self._lift(d)
        if d.is_zero():
            raise ZeroDivisionError("chia cho đa thức 0")
        r = list(self.coeffs)
        q = [self.F(0)] * max(0, len(r) - len(d.coeffs) + 1)
        inv_lead = d.coeffs[-1].inverse()
        while len(r) >= len(d.coeffs) and r:
            shift = len(r) - len(d.coeffs)
            c = r[-1] * inv_lead
            q[shift] = c
            for k, dc in enumerate(d.coeffs):
                r[shift + k] = r[shift + k] - c * dc
            while r and r[-1] == 0:
                r.pop()
        return Poly(self.F, q), Poly(self.F, r)

    def __floordiv__(self, d):
        return self.divmod(d)[0]

    def __mod__(self, d):
        return self.divmod(d)[1]

    def __call__(self, x):
        """Đánh giá tại x bằng lược đồ Horner."""
        acc = self.F(0) if not hasattr(x, "q") else x * 0
        for c in reversed(self.coeffs):
            acc = acc * x + c
        return acc

    def evaluate(self, x):
        return self(x)

    def scale(self, k):
        """Trả về f(k·X) (dùng cho coset)."""
        out, kp = [], self.F(1)
        for c in self.coeffs:
            out.append(c * kp)
            kp = kp * k
        return Poly(self.F, out)


def lagrange_interpolate(F, xs, ys):
    """Nội suy Lagrange: đa thức bậc < len(xs) đi qua các điểm (xs[i], ys[i])."""
    xs = [F(x) for x in xs]
    ys = [F(y) for y in ys]
    if len(set(int(x) for x in xs)) != len(xs):
        raise ValueError("các hoành độ phải phân biệt")
    result = Poly.zero(F)
    for i, xi in enumerate(xs):
        num, den = Poly.const(F, 1), F(1)
        for j, xj in enumerate(xs):
            if j != i:
                num = num * Poly(F, [-xj, 1])
                den = den * (xi - xj)
        result = result + num * (ys[i] / den)
    return result


def vanishing_poly(F, n):
    """Đa thức triệt tiêu Z_H(X) = X^n - 1 của nhóm con H = <ω> bậc n."""
    return Poly(F, [-1] + [0] * (n - 1) + [1])


def domain(F, n):
    """Miền H = [1, ω, ω^2, ..., ω^(n-1)] với ω = F.root_of_unity(n)."""
    w = F.root_of_unity(n)
    out, x = [], F(1)
    for _ in range(n):
        out.append(x)
        x = x * w
    return out


def lagrange_basis(F, n, i):
    """Đa thức cơ sở Lagrange L_i trên H: L_i(ω^j) = 1 nếu j = i, 0 nếu khác.

    Dạng đóng: L_i(X) = ω^i (X^n - 1) / (n (X - ω^i)).
    """
    H = domain(F, n)
    wi = H[i]
    q, r = vanishing_poly(F, n).divmod(Poly(F, [-wi, 1]))
    assert r.is_zero()
    return q * (wi / F(n))
