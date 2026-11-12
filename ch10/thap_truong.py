"""Tháp trường F_p -> F_{p^2} -> F_{p^6} -> F_{p^12} (cùng dạng với BLS12-381 trong arkworks):
    F_{p^2}  = F_p[u]/(u^2 + 1)             cần p ≡ 3 (mod 4)
    F_{p^6}  = F_{p^2}[v]/(v^3 - ξ)         ξ = 1 + u, cần ξ không là lập phương trong F_{p^2}
    F_{p^12} = F_{p^6}[w]/(w^2 - v)         cần v không là bình phương trong F_{p^6}
                                            (tương đương ξ không là bình phương trong F_{p^2})
Hệ quả: w^6 = ξ, tức F_{p^12} = F_{p^2}[w]/(w^6 - ξ).
(arkworks ark-bls12-381-0.6.0: fields/fq2.rs:13 NONRESIDUE = -1; fq6.rs:13-14
NONRESIDUE = U + 1; fq12.rs:13 NONRESIDUE = (0, 1, 0) = v.)

Phần tử: F_{p^2} là tuple (a0, a1) số nguyên mod p; F_{p^6} là tuple 3 phần tử F_{p^2};
F_{p^12} là tuple 2 phần tử F_{p^6}. DEM["mul"] đếm số phép nhân trong F_p.
Phép nhân dùng Karatsuba ở mọi tầng (3 / 6 / 3 phép nhân tầng dưới) -> 54 phép nhân F_p
cho một phép nhân F_{p^12}; bản "trường học" (karatsuba=False) cần 4·9·4 = 144.
"""

DEM = {"mul": 0, "inv": 0}


class Thap:
    def __init__(self, p, xi=(1, 1)):
        self.p = p
        self.xi = xi

    # ---------------- F_{p^2} ----------------
    def f2_add(self, a, b):
        return ((a[0] + b[0]) % self.p, (a[1] + b[1]) % self.p)

    def f2_sub(self, a, b):
        return ((a[0] - b[0]) % self.p, (a[1] - b[1]) % self.p)

    def f2_neg(self, a):
        return ((-a[0]) % self.p, (-a[1]) % self.p)

    def f2_mul(self, a, b, karatsuba=True):
        p = self.p
        if karatsuba:                       # 3 phép nhân F_p
            v0, v1 = a[0] * b[0], a[1] * b[1]
            m = (a[0] + a[1]) * (b[0] + b[1])
            DEM["mul"] += 3
            return ((v0 - v1) % p, (m - v0 - v1) % p)
        DEM["mul"] += 4
        return ((a[0] * b[0] - a[1] * b[1]) % p, (a[0] * b[1] + a[1] * b[0]) % p)

    def f2_mul_xi(self, a):
        """Nhân với ξ = 1 + u: (a0 + a1 u)(1 + u) = (a0 - a1) + (a0 + a1) u, không tốn phép nhân."""
        if self.xi == (1, 1):
            return ((a[0] - a[1]) % self.p, (a[0] + a[1]) % self.p)
        return self.f2_mul(a, self.xi)

    def f2_inv(self, a):
        """(a0 + a1 u)^{-1} = (a0 - a1 u)/(a0^2 + a1^2)."""
        n = (a[0] * a[0] + a[1] * a[1]) % self.p
        ni = pow(n, -1, self.p)
        DEM["mul"] += 4                     # a0^2, a1^2, a0·ni, a1·ni
        DEM["inv"] += 1                     # một nghịch đảo trong F_p
        return (a[0] * ni % self.p, (-a[1]) * ni % self.p)

    # ---------------- F_{p^6} ----------------
    def f6_add(self, a, b):
        return tuple(self.f2_add(x, y) for x, y in zip(a, b))

    def f6_sub(self, a, b):
        return tuple(self.f2_sub(x, y) for x, y in zip(a, b))

    def f6_mul(self, a, b, karatsuba=True):
        A, S, M, X = self.f2_add, self.f2_sub, lambda x, y: self.f2_mul(x, y, karatsuba), self.f2_mul_xi
        if karatsuba:                       # 6 phép nhân F_{p^2}
            v0, v1, v2 = M(a[0], b[0]), M(a[1], b[1]), M(a[2], b[2])
            c0 = A(v0, X(S(S(M(A(a[1], a[2]), A(b[1], b[2])), v1), v2)))
            c1 = A(S(S(M(A(a[0], a[1]), A(b[0], b[1])), v0), v1), X(v2))
            c2 = A(S(S(M(A(a[0], a[2]), A(b[0], b[2])), v0), v2), v1)
            return (c0, c1, c2)
        # trường học: 9 phép nhân F_{p^2}, rút gọn v^3 = ξ, v^4 = ξ v
        t = [[M(a[i], b[j]) for j in range(3)] for i in range(3)]
        c0 = A(t[0][0], X(A(t[1][2], t[2][1])))
        c1 = A(A(t[0][1], t[1][0]), X(t[2][2]))
        c2 = A(A(t[0][2], t[2][0]), t[1][1])
        return (c0, c1, c2)

    def f6_mul_v(self, a):
        """(a0 + a1 v + a2 v^2)·v = ξ a2 + a0 v + a1 v^2."""
        return (self.f2_mul_xi(a[2]), a[0], a[1])

    def f6_inv(self, a):
        M, S, X, A = self.f2_mul, self.f2_sub, self.f2_mul_xi, self.f2_add
        a0, a1, a2 = a
        c0 = S(M(a0, a0), X(M(a1, a2)))
        c1 = S(X(M(a2, a2)), M(a0, a1))
        c2 = S(M(a1, a1), M(a0, a2))
        F = A(M(a0, c0), X(A(M(a2, c1), M(a1, c2))))
        Fi = self.f2_inv(F)
        return (M(c0, Fi), M(c1, Fi), M(c2, Fi))

    # ---------------- F_{p^12} ----------------
    def f12_add(self, a, b):
        return (self.f6_add(a[0], b[0]), self.f6_add(a[1], b[1]))

    def f12_sub(self, a, b):
        return (self.f6_sub(a[0], b[0]), self.f6_sub(a[1], b[1]))

    def f12_mul(self, a, b, karatsuba=True):
        M = lambda x, y: self.f6_mul(x, y, karatsuba)
        if karatsuba:                       # 3 phép nhân F_{p^6}
            t0, t1 = M(a[0], b[0]), M(a[1], b[1])
            c1 = self.f6_sub(self.f6_sub(M(self.f6_add(a[0], a[1]), self.f6_add(b[0], b[1])), t0), t1)
            return (self.f6_add(t0, self.f6_mul_v(t1)), c1)
        t00, t01, t10, t11 = M(a[0], b[0]), M(a[0], b[1]), M(a[1], b[0]), M(a[1], b[1])
        return (self.f6_add(t00, self.f6_mul_v(t11)), self.f6_add(t01, t10))

    def f12_inv(self, a):
        """(a0 + a1 w)^{-1} = (a0 - a1 w)/(a0^2 - v a1^2)."""
        d = self.f6_sub(self.f6_mul(a[0], a[0]), self.f6_mul_v(self.f6_mul(a[1], a[1])))
        di = self.f6_inv(d)
        z = ((0, 0),) * 3
        return (self.f6_mul(a[0], di), self.f6_mul(self.f6_sub(z, a[1]), di))

    def f12_pow(self, a, e):
        r = self.one12()
        for bit in bin(e)[2:]:
            r = self.f12_mul(r, r)
            if bit == "1":
                r = self.f12_mul(r, a)
        return r

    # ---------------- tiện ích ----------------
    def zero12(self):
        return (((0, 0),) * 3,) * 2

    def one12(self):
        return (((1, 0), (0, 0), (0, 0)), ((0, 0),) * 3)

    def w(self):
        return (((0, 0),) * 3, ((1, 0), (0, 0), (0, 0)))

    def ngau_nhien12(self, rng):
        return tuple(tuple((rng.randrange(self.p), rng.randrange(self.p)) for _ in range(3))
                     for _ in range(2))

    def he_so(self, a):
        """12 hệ số F_p theo cơ sở (1, u)·(1, v, v^2)·(1, w)."""
        return [c for f6 in a for f2 in f6 for c in f2]


def f2_pow(p, a, e):
    T = Thap(p)
    r = (1, 0)
    for bit in bin(e)[2:]:
        r = T.f2_mul(r, r)
        if bit == "1":
            r = T.f2_mul(r, a)
    return r


def xi_hop_le(p, xi=(1, 1)):
    """(ξ không là bình phương, ξ không là lập phương) trong F_{p^2}, theo tiêu chuẩn Euler."""
    q = p * p - 1
    return f2_pow(p, xi, q // 2) != (1, 0), f2_pow(p, xi, q // 3) != (1, 0)
