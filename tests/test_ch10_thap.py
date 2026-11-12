"""Kiểm tháp F_p -> F_{p^2} -> F_{p^6} -> F_{p^12} của Ch.10 (ch10/thap_truong.py)."""
import os
import random
import sys

import _path  # noqa: F401

sys.path.insert(0, os.path.join(_path.ROOT, "ch10"))
from thap_truong import Thap, DEM, xi_hop_le  # noqa: E402


def _kiem_tien_de(T, rng, so_mau):
    M, A, S = T.f12_mul, T.f12_add, T.f12_sub
    one, zero = T.one12(), T.zero12()
    for _ in range(so_mau):
        x, y, z = (T.ngau_nhien12(rng) for _ in range(3))
        assert M(x, y) == M(y, x)
        assert M(M(x, y), z) == M(x, M(y, z))
        assert M(x, A(y, z)) == A(M(x, y), M(x, z))
        assert A(x, y) == A(y, x) and S(A(x, y), y) == x
        assert M(x, one) == x and A(x, zero) == x
        assert M(x, y) == T.f12_mul(x, y, karatsuba=False)
        if x != zero:
            assert M(x, T.f12_inv(x)) == one


def test_chon_p():
    assert xi_hop_le(19) == (True, True)
    assert xi_hop_le(7) == (False, True)          # ξ chính phương -> w^2 - v khả quy
    assert xi_hop_le(43) == (True, False)         # ξ lập phương -> v^3 - ξ khả quy


def test_tien_de_19_va_67():
    rng = random.Random(449)
    for p in (19, 67):
        _kiem_tien_de(Thap(p), rng, 60)


def test_cau_truc():
    T = Thap(19)
    w = T.w()
    assert T.f12_pow(w, 6) == (((1, 1), (0, 0), (0, 0)), ((0, 0),) * 3)
    assert all(T.f12_pow(w, 19 ** k) != w for k in range(1, 12))
    assert T.f12_pow(w, 19 ** 12) == w


def test_dem_phep_nhan():
    T = Thap(19)
    rng = random.Random(1)
    a, b = T.ngau_nhien12(rng), T.ngau_nhien12(rng)
    DEM["mul"] = 0
    T.f12_mul(a, b)
    assert DEM["mul"] == 54
    DEM["mul"] = 0
    T.f12_mul(a, b, karatsuba=False)
    assert DEM["mul"] == 144


if __name__ == "__main__":
    _path.run_module(globals())
