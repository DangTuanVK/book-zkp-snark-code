"""Trường tháp nhị phân F_2 ⊂ F_4 ⊂ F_16 ⊂ F_256 theo cấu trúc Wiedemann (Binius).

    T_0 = F_2,
    T_{ι+1} = T_ι[X_ι] / (X_ι^2 + X_{ι-1}·X_ι + 1),   quy ước X_{-1} = 1.
    |T_ι| = 2^{2^ι}:  T_1 = F_4, T_2 = F_16, T_3 = F_256.

Biểu diễn: phần tử của T_ι là số nguyên 2^ι bit. Phần tử của T_{ι+1} là
a = a_0 + a_1·X_ι với a_0 = nửa bit THẤP, a_1 = nửa bit CAO (đều thuộc T_ι).
Hệ quả: bit thứ i ứng với đơn thức Π_{j : bit j của i bằng 1} X_j (cơ sở đa
tuyến tính), ví dụ trong F_16: bit0 = 1, bit1 = X_0, bit2 = X_1, bit3 = X_0·X_1.

Cộng = XOR. Nhân (c = X_{ι-1} ∈ T_ι, dùng X_ι^2 = c·X_ι + 1):
    (a_0 + a_1 X)(b_0 + b_1 X) = (a_0 b_0 + a_1 b_1) + (a_0 b_1 + a_1 b_0 + a_1 b_1·c)·X
Nghịch đảo: δ = a_0^2 + a_0 a_1 c + a_1^2 ∈ T_ι,
    a^{-1} = ((a_0 + a_1 c) + a_1·X)·δ^{-1}.

Nguồn: Diamond–Posen, ePrint 2023/1784 (cấu trúc tháp của Wiedemann 1988);
công thức cần đối chiếu khi có bản PDF.
"""

MAX_LEVEL = 3


def bits_of(lvl):
    return 1 << lvl


def _check(a, lvl):
    if not 0 <= lvl <= MAX_LEVEL + 3:
        raise ValueError("mức tháp không hợp lệ")
    if not 0 <= a < (1 << bits_of(lvl)):
        raise ValueError("%d không thuộc T_%d" % (a, lvl))


def gen(lvl):
    """X_{lvl-1} — phần tử sinh 'mới' của T_lvl (lvl ≥ 1); X_{-1} = 1 ở T_0."""
    return 1 if lvl == 0 else 1 << (bits_of(lvl) // 2)


def add(a, b):
    return a ^ b


def mul(a, b, lvl):
    """Nhân trong T_lvl, đệ quy kiểu Karatsuba (3 phép nhân ở tầng dưới)."""
    if lvl == 0:
        return a & b
    h = bits_of(lvl) // 2
    mask = (1 << h) - 1
    a0, a1, b0, b1 = a & mask, a >> h, b & mask, b >> h
    p00 = mul(a0, b0, lvl - 1)
    p11 = mul(a1, b1, lvl - 1)
    pm = mul(a0 ^ a1, b0 ^ b1, lvl - 1) ^ p00 ^ p11      # = a0 b1 + a1 b0
    c = gen(lvl - 1)
    lo = p00 ^ p11
    hi = pm ^ mul(p11, c, lvl - 1)
    return lo | (hi << h)


def square(a, lvl):
    return mul(a, a, lvl)


def power(a, e, lvl):
    r = 1
    while e:
        if e & 1:
            r = mul(r, a, lvl)
        a = mul(a, a, lvl)
        e >>= 1
    return r


def inv(a, lvl):
    """Nghịch đảo đệ quy qua chuẩn δ = a0^2 + a0 a1 c + a1^2."""
    if a == 0:
        raise ZeroDivisionError("0 không có nghịch đảo")
    if lvl == 0:
        return 1
    h = bits_of(lvl) // 2
    mask = (1 << h) - 1
    a0, a1 = a & mask, a >> h
    c = gen(lvl - 1)
    a1c = mul(a1, c, lvl - 1)
    delta = mul(a0, a0, lvl - 1) ^ mul(a0, a1c, lvl - 1) ^ mul(a1, a1, lvl - 1)
    di = inv(delta, lvl - 1)
    return mul(a0 ^ a1c, di, lvl - 1) | (mul(a1, di, lvl - 1) << h)


def fmt(a, lvl):
    """Viết phần tử dưới dạng tổng các đơn thức X_j (cơ sở đa tuyến tính)."""
    if a == 0:
        return "0"
    terms = []
    for i in range(bits_of(lvl)):
        if a >> i & 1:
            mono = [("X%d" % j) for j in range(lvl) if i >> j & 1]
            terms.append("·".join(mono) if mono else "1")
    return " + ".join(terms)


def mul_trace(a, b, lvl):
    """Các bước nhân 'trường học' một tầng (4 tích con) — dùng cho ví dụ tính tay."""
    h = bits_of(lvl) // 2
    mask = (1 << h) - 1
    a0, a1, b0, b1 = a & mask, a >> h, b & mask, b >> h
    c = gen(lvl - 1)
    p00 = mul(a0, b0, lvl - 1)
    p11 = mul(a1, b1, lvl - 1)
    p01 = mul(a0, b1, lvl - 1) ^ mul(a1, b0, lvl - 1)
    p11c = mul(p11, c, lvl - 1)
    lo, hi = p00 ^ p11, p01 ^ p11c
    return {"a0": a0, "a1": a1, "b0": b0, "b1": b1, "c": c,
            "a0b0": p00, "a1b1": p11, "a0b1+a1b0": p01, "a1b1*c": p11c,
            "lo": lo, "hi": hi, "ket_qua": lo | (hi << h)}


class TowerElement:
    """Phần tử của T_lvl với toán tử + (XOR), *, inverse, ** cho mã dễ đọc."""

    __slots__ = ("v", "lvl")

    def __init__(self, v, lvl):
        _check(v, lvl)
        self.v, self.lvl = v, lvl

    def _o(self, o):
        if isinstance(o, TowerElement):
            if o.lvl != self.lvl:
                raise ValueError("khác mức tháp (dùng .lift())")
            return o.v
        return TowerElement(o, self.lvl).v

    def __add__(self, o):
        return TowerElement(self.v ^ self._o(o), self.lvl)

    __sub__ = __radd__ = __rsub__ = __add__

    def __mul__(self, o):
        return TowerElement(mul(self.v, self._o(o), self.lvl), self.lvl)

    __rmul__ = __mul__

    def __pow__(self, e):
        if e < 0:
            return self.inverse() ** (-e)
        return TowerElement(power(self.v, e, self.lvl), self.lvl)

    def inverse(self):
        return TowerElement(inv(self.v, self.lvl), self.lvl)

    def __truediv__(self, o):
        return self * TowerElement(self._o(o), self.lvl).inverse()

    def lift(self, lvl):
        """Nhúng vào tầng cao hơn (giữ nguyên chuỗi bit — tháp lồng nhau)."""
        if lvl < self.lvl:
            raise ValueError("chỉ nâng lên tầng cao hơn")
        return TowerElement(self.v, lvl)

    def __eq__(self, o):
        return isinstance(o, TowerElement) and (self.v, self.lvl) == (o.v, o.lvl)

    def __hash__(self):
        return hash((self.v, self.lvl))

    def __repr__(self):
        return "0b{:0{w}b}".format(self.v, w=bits_of(self.lvl))


def T(lvl):
    """Hàm tạo phần tử của T_lvl: F16 = T(2); F16(0b0110)."""
    return lambda v: TowerElement(v, lvl)
