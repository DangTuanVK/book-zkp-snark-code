"""Kiểm tham số BN254 / BLS12-381 suy từ u khớp arkworks 0.6.0 (Ch.10)."""
import os
import sys

import _path  # noqa: F401

sys.path.insert(0, os.path.join(_path.ROOT, "ch10"))
import ho_duong_cong as H  # noqa: E402


def test_bn254():
    p, r, t = H.bn(H.U_BN254)
    assert (p, r) == (H.ARK_P_BN254, H.ARK_R_BN254)
    assert H.miller_rabin(p) and H.miller_rabin(r)
    assert H.bac_nhung(p, r) == 12 and H.two_adicity(r - 1) == 28
    assert p + 1 - t == r                               # đồng hệ số 1
    assert H.naf_gia_tri(H.ARK_ATE_LOOP_BN254) == 6 * H.U_BN254 + 2


def test_bls12_381():
    u = H.U_BLS12_381
    p, r, t = H.bls12(u)
    assert (p, r) == (H.ARK_P_BLS12_381, H.ARK_R_BLS12_381)
    assert H.miller_rabin(p) and H.miller_rabin(r)
    assert H.bac_nhung(p, r) == 12 and H.two_adicity(r - 1) == 32
    assert (p + 1 - t) % r == 0 and (p + 1 - t) // r == H.ARK_H1_BLS12_381
    assert (p ** 12 - 1) // r == (p ** 6 - 1) * (p ** 2 + 1) * ((p ** 4 - p ** 2 + 1) // r)


def test_miller_rabin_nho():
    assert [n for n in range(2, 60) if H.miller_rabin(n)] == \
        [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59]
    assert not H.miller_rabin(3215031751)               # hợp số Carmichael mạnh với cơ sở 2,3,5,7
    assert not H.miller_rabin(391)                      # 17·23


if __name__ == "__main__":
    _path.run_module(globals())
