"""Kiểm snarklib/montgomery.py: REDC, nhân Montgomery, CIOS (Ch.10)."""
import os
import random
import sys

import _path  # noqa: F401
from snarklib.montgomery import Montgomery

sys.path.insert(0, os.path.join(_path.ROOT, "ch10"))
from ho_duong_cong import ARK_P_BLS12_381, ARK_P_BN254  # noqa: E402


def test_449_vet_can():
    M = Montgomery(449)
    assert M.R == 512 and (M.p * M.p_prime) % M.R == M.R - 1
    for a in range(449):
        am = M.to_mont(a)
        assert am == a * M.R % 449 and M.from_mont(am) == a
        for b in range(0, 449, 7):
            assert M.from_mont(M.mul(am, M.to_mont(b))) == a * b % 449


def test_bls12_381_va_bn254():
    rng = random.Random(449)
    for p, s in ((ARK_P_BLS12_381, 6), (ARK_P_BN254, 4)):
        M = Montgomery(p, word=64)
        assert M.s == s and M.R == 1 << (64 * s)
        for _ in range(300):
            x, y = rng.randrange(p), rng.randrange(p)
            xm, ym = M.to_mont(x), M.to_mont(y)
            z, nhan = M.mul_cios(xm, ym)
            assert z == M.mul(xm, ym) == x * y * M.R % p
            assert nhan == 2 * s * s + s
        e = rng.randrange(p)
        assert M.from_mont(M.pow(M.to_mont(3), e)) == pow(3, e, p)


def test_bien_redc():
    M = Montgomery(449)
    T = M.p * M.R - 1                      # giá trị lớn nhất cho phép
    assert M.redc(T) == T * pow(M.R, -1, 449) % 449


if __name__ == "__main__":
    _path.run_module(globals())
