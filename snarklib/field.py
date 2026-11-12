"""Trường hữu hạn F_p và trường mở rộng bậc 2 F_q[i]/(i^2 + 1).

Cách dùng:
    F = GF(449)          # lớp phần tử của F_449
    a = F(3); b = a**28  # b = 349
    F.generator()        # phần tử sinh nhỏ nhất của F_p^*
    F.root_of_unity(16)  # g^((p-1)/16)

    K = GF2(5387)        # F_{q^2} = F_q[i]/(i^2+1), cần q ≡ 3 (mod 4)
    i = K(0, 1); i*i == K(-1)
"""

_CACHE = {}


def is_prime(n):
    """Kiểm nguyên tố bằng chia thử (đủ cho các số nhỏ trong sách)."""
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def factorize(n):
    """Phân tích n thành thừa số nguyên tố: trả về dict {p: số mũ}."""
    out = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


class FieldElement:
    """Phần tử của F_p. Lớp con do GF(p) tạo ra mang thuộc tính p."""

    p = None
    __slots__ = ("v",)

    def __init__(self, v):
        if isinstance(v, FieldElement):
            if v.p != self.p:
                raise ValueError("khác trường")
            v = v.v
        self.v = v % self.p

    # --- chuyển kiểu ---
    def _coerce(self, other):
        if isinstance(other, FieldElement):
            if other.p != self.p:
                raise ValueError("phép toán giữa hai trường khác nhau")
            return other.v
        if isinstance(other, int):
            return other % self.p
        return NotImplemented

    def __int__(self):
        return self.v

    def __index__(self):
        return self.v

    def __repr__(self):
        return str(self.v)

    def __hash__(self):
        return hash((self.p, self.v))

    # --- số học ---
    def __add__(self, o):
        o = self._coerce(o)
        if o is NotImplemented:
            return o
        return type(self)(self.v + o)

    __radd__ = __add__

    def __sub__(self, o):
        o = self._coerce(o)
        if o is NotImplemented:
            return o
        return type(self)(self.v - o)

    def __rsub__(self, o):
        o = self._coerce(o)
        if o is NotImplemented:
            return o
        return type(self)(o - self.v)

    def __neg__(self):
        return type(self)(-self.v)

    def __mul__(self, o):
        o = self._coerce(o)
        if o is NotImplemented:
            return o
        return type(self)(self.v * o)

    __rmul__ = __mul__

    def inverse(self):
        """Nghịch đảo nhân theo định lý Fermat nhỏ: a^(p-2)."""
        if self.v == 0:
            raise ZeroDivisionError("0 không có nghịch đảo")
        return type(self)(pow(self.v, self.p - 2, self.p))

    def __truediv__(self, o):
        o = self._coerce(o)
        if o is NotImplemented:
            return o
        return self * type(self)(o).inverse()

    def __rtruediv__(self, o):
        o = self._coerce(o)
        if o is NotImplemented:
            return o
        return type(self)(o) * self.inverse()

    def __pow__(self, e):
        if e < 0:
            return self.inverse() ** (-e)
        return type(self)(pow(self.v, e, self.p))

    # --- so sánh ---
    def __eq__(self, o):
        o = self._coerce(o)
        if o is NotImplemented:
            return False
        return self.v == o

    def __ne__(self, o):
        return not self.__eq__(o)

    def __bool__(self):
        return self.v != 0

    # --- thặng dư bậc hai, căn bậc hai ---
    def legendre(self):
        """Ký hiệu Legendre: 1 (thặng dư khác 0), -1 (không thặng dư), 0."""
        if self.v == 0:
            return 0
        t = pow(self.v, (self.p - 1) // 2, self.p)
        return 1 if t == 1 else -1

    def is_qr(self):
        """True nếu là thặng dư bậc hai (kể cả 0)."""
        return self.legendre() >= 0

    def sqrt(self):
        """Căn bậc hai bằng Tonelli–Shanks; trả về căn NHỎ HƠN (theo số nguyên).

        Ném ValueError nếu không phải thặng dư bậc hai.
        """
        p = self.p
        if self.v == 0:
            return type(self)(0)
        if self.legendre() != 1:
            raise ValueError("%d không là thặng dư bậc hai mod %d" % (self.v, p))
        # p - 1 = Q * 2^S với Q lẻ
        Q, S = p - 1, 0
        while Q % 2 == 0:
            Q //= 2
            S += 1
        # z: một số không thặng dư bất kỳ
        z = 2
        while pow(z, (p - 1) // 2, p) != p - 1:
            z += 1
        M, c = S, pow(z, Q, p)
        t, R = pow(self.v, Q, p), pow(self.v, (Q + 1) // 2, p)
        while t != 1:
            # tìm i nhỏ nhất với t^(2^i) = 1
            i, t2 = 0, t
            while t2 != 1:
                t2 = t2 * t2 % p
                i += 1
            b = pow(c, 1 << (M - i - 1), p)
            M, c = i, b * b % p
            t, R = t * c % p, R * b % p
        return type(self)(min(R, p - R))

    def order(self):
        """Bậc nhân của phần tử trong F_p^*."""
        if self.v == 0:
            raise ValueError("0 không thuộc F_p^*")
        n = self.p - 1
        for q, e in factorize(self.p - 1).items():
            for _ in range(e):
                if pow(self.v, n // q, self.p) == 1:
                    n //= q
                else:
                    break
        return n

    # --- hàm mức trường (gọi qua lớp: F.generator()) ---
    @classmethod
    def zero(cls):
        return cls(0)

    @classmethod
    def one(cls):
        return cls(1)

    @classmethod
    def generator(cls):
        """Phần tử sinh nhỏ nhất của nhóm nhân F_p^*."""
        p = cls.p
        primes = list(factorize(p - 1))
        for g in range(2, p):
            if all(pow(g, (p - 1) // q, p) != 1 for q in primes):
                return cls(g)
        return cls(1)  # p = 2

    @classmethod
    def two_adicity(cls):
        """Số mũ s lớn nhất với 2^s | p - 1."""
        n, s = cls.p - 1, 0
        while n % 2 == 0:
            n //= 2
            s += 1
        return s

    @classmethod
    def root_of_unity(cls, n):
        """Căn nguyên thủy bậc n của đơn vị: g^((p-1)/n), g là phần tử sinh nhỏ nhất."""
        if (cls.p - 1) % n != 0:
            raise ValueError("n = %d không chia hết p - 1 = %d" % (n, cls.p - 1))
        w = cls.generator() ** ((cls.p - 1) // n)
        assert w.order() == n
        return w

    @classmethod
    def elements(cls):
        """Duyệt mọi phần tử 0, 1, ..., p-1."""
        for v in range(cls.p):
            yield cls(v)


def GF(p):
    """Tạo (và lưu đệm) lớp phần tử của trường F_p, p nguyên tố."""
    if p in _CACHE:
        return _CACHE[p]
    if not is_prime(p):
        raise ValueError("%d không nguyên tố" % p)
    cls = type("F%d" % p, (FieldElement,), {"p": p, "__slots__": ()})
    _CACHE[p] = cls
    return cls


# ---------------------------------------------------------------------------
# Trường mở rộng bậc 2: F_{q^2} = F_q[i]/(i^2 + 1), hợp lệ khi q ≡ 3 (mod 4)
# (khi đó -1 không là thặng dư bậc hai nên i^2 + 1 bất khả quy).
# ---------------------------------------------------------------------------

class Fq2Element:
    """Phần tử a + b·i của F_q[i]/(i^2+1). Lớp con do GF2(q) tạo ra."""

    q = None
    __slots__ = ("a", "b")

    def __init__(self, a, b=0):
        q = self.q
        self.a = int(a) % q
        self.b = int(b) % q

    def _coerce(self, o):
        if isinstance(o, Fq2Element):
            if o.q != self.q:
                raise ValueError("khác trường")
            return o
        if isinstance(o, (int, FieldElement)):
            return type(self)(int(o), 0)
        return NotImplemented

    def __repr__(self):
        if self.b == 0:
            return "%d" % self.a
        return "%d + %d·i" % (self.a, self.b)

    def __hash__(self):
        return hash((self.q, self.a, self.b))

    def __eq__(self, o):
        o = self._coerce(o)
        if o is NotImplemented:
            return False
        return self.a == o.a and self.b == o.b

    def __ne__(self, o):
        return not self.__eq__(o)

    def __bool__(self):
        return self.a != 0 or self.b != 0

    def __add__(self, o):
        o = self._coerce(o)
        if o is NotImplemented:
            return o
        return type(self)(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __sub__(self, o):
        o = self._coerce(o)
        if o is NotImplemented:
            return o
        return type(self)(self.a - o.a, self.b - o.b)

    def __rsub__(self, o):
        return self._coerce(o) - self

    def __neg__(self):
        return type(self)(-self.a, -self.b)

    def __mul__(self, o):
        o = self._coerce(o)
        if o is NotImplemented:
            return o
        # (a + bi)(c + di) = (ac - bd) + (ad + bc) i
        return type(self)(self.a * o.a - self.b * o.b, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    def conjugate(self):
        """Liên hợp a - b·i (cũng là ánh xạ Frobenius x -> x^q)."""
        return type(self)(self.a, -self.b)

    def norm(self):
        """Chuẩn a^2 + b^2 thuộc F_q."""
        return (self.a * self.a + self.b * self.b) % self.q

    def inverse(self):
        n = self.norm()
        if n == 0:
            raise ZeroDivisionError("0 không có nghịch đảo")
        ninv = pow(n, self.q - 2, self.q)
        return type(self)(self.a * ninv, -self.b * ninv)

    def __truediv__(self, o):
        o = self._coerce(o)
        if o is NotImplemented:
            return o
        return self * o.inverse()

    def __pow__(self, e):
        if e < 0:
            return self.inverse() ** (-e)
        result, base = type(self)(1), self
        while e:
            if e & 1:
                result = result * base
            base = base * base
            e >>= 1
        return result


_CACHE2 = {}


def GF2(q):
    """Tạo lớp phần tử của F_{q^2} = F_q[i]/(i^2+1); yêu cầu q nguyên tố, q ≡ 3 mod 4."""
    if q in _CACHE2:
        return _CACHE2[q]
    if not is_prime(q) or q % 4 != 3:
        raise ValueError("cần q nguyên tố và q ≡ 3 (mod 4)")
    cls = type("F%d_2" % q, (Fq2Element,), {"q": q, "__slots__": ()})
    _CACHE2[q] = cls
    return cls
