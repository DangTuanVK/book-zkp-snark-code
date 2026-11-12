"""Kiểm r1cs.py trên ví dụ a·b = N = 391."""
import _path  # noqa: F401
from snarklib.field import GF
from snarklib.r1cs import R1CS, factor_r1cs

F = GF(449)


def test_r1cs():
    r = factor_r1cs(F)
    assert r.is_satisfied([391], [17, 23])
    assert r.is_satisfied([391], [28, 30])      # nghiệm giả: 840 ≡ 391 mod 449
    assert not r.is_satisfied([391], [17, 22])
    # R1CS 2 ràng buộc: v = a·a, v·a = c (c = a^3)
    A = [[0, 0, 1, 0], [0, 0, 0, 1]]
    B = [[0, 0, 1, 0], [0, 0, 1, 0]]
    C = [[0, 0, 0, 1], [0, 1, 0, 0]]
    cube = R1CS(F, A, B, C, n_public=1)          # z = (1, c, a, v)
    assert cube.is_satisfied([27], [3, 9]) and not cube.is_satisfied([28], [3, 9])


if __name__ == "__main__":
    _path.run_module(globals())
