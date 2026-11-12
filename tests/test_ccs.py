"""Kiểm ccs.py: CCS từ R1CS và CCS bậc 3."""
import _path  # noqa: F401
from snarklib.field import GF
from snarklib.r1cs import factor_r1cs
from snarklib.ccs import CCS, from_r1cs

F = GF(449)


def test_ccs():
    r = factor_r1cs(F)
    s = from_r1cs(r)
    assert s.degree() == 2
    for w in ([17, 23], [28, 30], [17, 22], [1, 1]):
        assert s.is_satisfied([391], w) == r.is_satisfied([391], w)
    # CCS bậc 3: a·b·c = N với z = (1, N, a, b, c)
    Ma = [[0, 0, 1, 0, 0]]
    Mb = [[0, 0, 0, 1, 0]]
    Mc = [[0, 0, 0, 0, 1]]
    MN = [[0, 1, 0, 0, 0]]
    s3 = CCS(F, [Ma, Mb, Mc, MN], [[0, 1, 2], [3]], [1, -1], n_public=1)
    assert s3.degree() == 3
    assert s3.is_satisfied([30], [2, 3, 5]) and not s3.is_satisfied([31], [2, 3, 5])


if __name__ == "__main__":
    _path.run_module(globals())
