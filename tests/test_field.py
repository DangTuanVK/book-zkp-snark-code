"""Kiểm field.py: F_p, Tonelli–Shanks, phần tử sinh, căn đơn vị, F_{q^2}."""
import _path  # noqa: F401
from snarklib.field import GF, GF2, is_prime, factorize


def test_arith():
    F = GF(449)
    a, b = F(17), F(23)
    assert a * b == 391
    assert F(28) * F(30) == 391           # nghiệm giả do tràn trường
    assert (a / b) * b == a
    assert a - b == 443 and -a == 432
    assert a ** 448 == 1 and a ** -1 * a == 1
    assert 5 - a == F(-12) and 3 / F(3) == 1
    for v in range(1, 449):
        assert F(v) * F(v).inverse() == 1


def test_sqrt_all():
    for p in (449, 5387, 17, 97, 257):
        F = GF(p)
        n_qr = 0
        for v in range(1, p):
            x = F(v)
            if x.is_qr():
                n_qr += 1
                s = x.sqrt()
                assert s * s == x and int(s) <= p - int(s)
            else:
                try:
                    x.sqrt()
                    assert False
                except ValueError:
                    pass
        assert n_qr == (p - 1) // 2


def test_generator_roots():
    F = GF(449)
    assert is_prime(449) and factorize(448) == {2: 6, 7: 1}
    assert F.generator() == 3 and F(3).order() == 448
    assert F.two_adicity() == 6
    assert F.root_of_unity(16) == 349 == F(3) ** 28
    for n in (2, 4, 7, 8, 14, 16, 32, 64, 448):
        w = F.root_of_unity(n)
        assert w ** n == 1 and w.order() == n


def test_fq2():
    q = 5387
    K = GF2(q)
    i = K(0, 1)
    assert i * i == K(-1)
    x = K(1234, 4321)
    assert x * x.inverse() == 1
    assert x ** (q * q - 1) == 1
    assert x ** q == x.conjugate()        # Frobenius
    try:
        GF2(449)                          # 449 ≡ 1 mod 4: không hợp lệ
        assert False
    except ValueError:
        pass


if __name__ == "__main__":
    _path.run_module(globals())
