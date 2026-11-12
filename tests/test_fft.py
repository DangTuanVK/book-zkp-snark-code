"""Kiểm fft.py: NTT/INTT, coset NTT khớp đánh giá trực tiếp."""
import random
import _path  # noqa: F401
from snarklib.field import GF
from snarklib.fft import ntt, intt, coset_ntt, coset_intt, naive_eval
from snarklib.poly import domain

F = GF(449)


def test_ntt_matches_naive():
    rng = random.Random(449)
    for n in (1, 2, 4, 8, 16, 32, 64):
        a = [F(rng.randrange(449)) for _ in range(n)]
        H = domain(F, n) if n > 1 else [F(1)]
        assert ntt(F, a) == naive_eval(F, a, H)
        assert intt(F, ntt(F, a)) == a


def test_coset():
    rng = random.Random(7)
    for n in (8, 16):
        a = [F(rng.randrange(449)) for _ in range(n)]
        k = F(5)
        pts = [k * h for h in domain(F, n)]
        assert coset_ntt(F, a, k) == naive_eval(F, a, pts)
        assert coset_intt(F, coset_ntt(F, a, k), k) == a


def test_other_field():
    G = GF(257)                           # 257 - 1 = 2^8
    a = list(range(1, 129))
    a = [G(x) for x in a]
    assert ntt(G, a) == naive_eval(G, a, domain(G, 128))
    try:
        ntt(F, [F(0)] * 128)              # F_449 không có căn bậc 128
        assert False
    except ValueError:
        pass


if __name__ == "__main__":
    _path.run_module(globals())
