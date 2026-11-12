"""Kiểm snarklib/pippenger.py: Pippenger = ngây thơ = tính trực tiếp, và số đếm (Ch.10)."""
import random

import _path  # noqa: F401
from snarklib.curve import toy_generator, R_TOY
from snarklib.pippenger import (NhomDiem, msm_pippenger, msm_ngay_tho, dem_pippenger,
                                dem_ngay_tho, chi_phi_ly_thuyet)


def test_pippenger_bang_ngay_tho():
    rng = random.Random(12)
    E, P, G = toy_generator()
    for n in (1, 2, 5, 16, 33):
        k = [rng.randrange(1, R_TOY) for _ in range(n)]
        pts = [ki * G for ki in k]
        for b in (9, 64):
            s = [rng.getrandbits(b) for _ in range(n)]
            s[0] = 0                                   # vô hướng 0 phải được xử lý
            mong_doi = (sum(a * c for a, c in zip(s, k)) % R_TOY) * G
            assert msm_ngay_tho(NhomDiem(), pts, s) == mong_doi
            for c in (1, 2, 3, 4, 7, 9):
                assert msm_pippenger(NhomDiem(), pts, s, c, bits=b) == mong_doi


def test_dem_khop_tinh_that():
    """Số đếm ở chế độ chỉ-đếm bằng số đếm khi tính trên điểm thật."""
    rng = random.Random(5)
    E, P, G = toy_generator()
    pts = [rng.randrange(1, R_TOY) * G for _ in range(20)]
    s = [rng.getrandbits(30) for _ in range(20)]
    for c in (2, 5):
        Gt = NhomDiem()
        msm_pippenger(Gt, pts, s, c, bits=30)
        d = dem_pippenger(s, c, 30)
        assert (Gt.dem.add, Gt.dem.dbl) == (d.add, d.dbl)
    Gt = NhomDiem()
    msm_ngay_tho(Gt, pts, s)
    d = dem_ngay_tho(s)
    assert (Gt.dem.add, Gt.dem.dbl) == (d.add, d.dbl)


def test_can_tren_ly_thuyet():
    rng = random.Random(7)
    for n in (16, 256):
        s = [rng.getrandbits(255) for _ in range(n)]
        for c in range(1, 11):
            d = dem_pippenger(s, c, 255)
            a, db = chi_phi_ly_thuyet(n, 255, c)
            assert d.add <= a and d.dbl == db


if __name__ == "__main__":
    _path.run_module(globals())
