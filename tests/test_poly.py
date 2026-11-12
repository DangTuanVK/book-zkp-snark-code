"""Kiểm poly.py: số học, chia có dư, nội suy, Z_H, L_i."""
import _path  # noqa: F401
from snarklib.field import GF
from snarklib.poly import Poly, lagrange_interpolate, vanishing_poly, domain, lagrange_basis

F = GF(449)


def test_arith_divmod():
    f = Poly(F, [5, 0, 3, 1])
    g = Poly(F, [-2, 1])
    q, r = f.divmod(g)
    assert q * g + r == f and r.degree() < g.degree()
    assert r == f(2)                      # định lý dư
    assert (f * g).degree() == 4 and (f - f).is_zero()
    assert (g ** 3)(F(5)) == F(27)


def test_interpolate():
    xs, ys = [1, 2, 3, 4], [17, 23, 391, 0]
    p = lagrange_interpolate(F, xs, ys)
    assert p.degree() <= 3
    assert all(p(x) == y for x, y in zip(xs, ys))


def test_vanishing_lagrange():
    for n in (4, 8, 16):
        H = domain(F, n)
        Z = vanishing_poly(F, n)
        assert all(Z(h) == 0 for h in H) and Z(F(2)) != 0
        total = Poly.zero(F)
        for i in range(n):
            L = lagrange_basis(F, n, i)
            assert L.degree() == n - 1
            assert all(L(H[j]) == (1 if i == j else 0) for j in range(n))
            total = total + L
        assert total == Poly.const(F, 1)  # Σ L_i = 1


if __name__ == "__main__":
    _path.run_module(globals())
