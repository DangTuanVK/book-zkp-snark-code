"""Test phần Python của Ch.11: mô hình bit-chính-xác của montmul.sv/ntt_bfly.sv và
bộ sinh vector (không cần VCS; mô phỏng thật chạy bằng rtl/run_all.sh)."""
import random
import _path  # noqa: F401
from snarklib.montgomery import Montgomery
from snarklib.field import GF
from snarklib import fft
import importlib.util, os

_spec = importlib.util.spec_from_file_location(
    "sinh_vector_rtl", os.path.join(_path.ROOT, "ch11", "sinh_vector_rtl.py"))
sv = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sv)


def redc_rtl(a, b, w, p, pinv):
    R = 1 << w
    T = a * b
    m = (T % R) * pinv % R
    U = T + m * p
    assert U % R == 0
    t = U >> w
    assert t < 2 * p
    return t - p if t >= p else t


def test_vi_du_tay_w16():
    M = Montgomery(449, word=16)
    assert M.p_prime == 0xF1BF and M.R_mod_p == 431
    a, b = M.to_mont(17), M.to_mont(23)
    assert (a, b) == (143, 35)
    assert redc_rtl(a, b, 16, 449, M.p_prime) == 146 == M.mul(a, b)
    assert M.from_mont(146) == 391


def test_cac_cau_hinh_ngau_nhien():
    rng = random.Random(7)
    for _, w, p, _ in sv.CAU_HINH:
        M = Montgomery(p, word=w)
        assert M.k == w and p < (1 << w) and p % 2 == 1
        for _ in range(300):
            a, b = rng.randrange(p), rng.randrange(p)
            y = redc_rtl(a, b, w, p, M.p_prime)
            assert y == M.mul(a, b)
            assert M.from_mont(y) == M.from_mont(a) * M.from_mont(b) % p


def test_bfly_vet_ntt_khop_fft():
    F = GF(449)
    rng = random.Random(3)
    M, vs = sv.bfly_vectors(449, 16, rng, (8, 16, 32), None, 50, F)
    assert len(vs) == 12 + 32 + 80 + 50 + 4   # vết NTT n=8,16,32 + ngẫu nhiên + biên
    for u, v, wm, x, y in vs:
        om = wm * pow(M.R, -1, 449) % 449      # twiddle dạng thường
        assert x == (u + v * om) % 449 and y == (u - v * om) % 449


def test_goldilocks_can_don_vi():
    GL = 2**64 - 2**32 + 1
    om = pow(7, (GL - 1) // 64, GL)
    assert pow(om, 64, GL) == 1 and pow(om, 32, GL) != 1


def test_so_lan_tru_cuoi_w16():
    M = Montgomery(449, word=16)
    dem = sum(1 for x in range(449) for z in range(449)
              if ((x * z + ((x * z * M.p_prime) & M.mask) * 449) >> 16) >= 449)
    assert dem == 155


if __name__ == "__main__":
    _path.run_module(globals())
