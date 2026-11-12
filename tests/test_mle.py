"""Kiểm mle.py: MLE trùng bảng trên hypercube, eq, fold."""
import random
import _path  # noqa: F401
from snarklib.field import GF
from snarklib.mle import eq, eq_table, evaluate, fold, evaluate_by_folding, bits

F = GF(449)


def test_hypercube():
    T = [F(v) for v in [17, 23, 391, 5, 0, 1, 2, 448]]
    for i in range(8):
        assert evaluate(F, T, bits(i, 3)) == T[i]
    assert bits(1, 2) == [0, 1] and bits(2, 2) == [1, 0]   # x_1 là bit cao


def test_eq():
    for i in range(8):
        for j in range(8):
            assert eq(F, bits(i, 3), bits(j, 3)) == (1 if i == j else 0)
    r = [F(3), F(100), F(7)]
    tab = eq_table(F, r)
    assert all(tab[i] == eq(F, bits(i, 3), r) for i in range(8))
    assert sum(tab, F(0)) == 1


def test_fold():
    rng = random.Random(1)
    for l in (1, 2, 3, 4, 5):
        T = [F(rng.randrange(449)) for _ in range(2 ** l)]
        r = [F(rng.randrange(449)) for _ in range(l)]
        assert evaluate(F, T, r) == evaluate_by_folding(F, T, r)
        # fold x_1 = r_1 rồi đánh giá phần còn lại = đánh giá toàn bộ
        if l > 1:
            assert evaluate(F, fold(F, T, r[0]), r[1:]) == evaluate(F, T, r)
    # đa tuyến tính: tuyến tính theo từng biến
    T = [F(rng.randrange(449)) for _ in range(8)]
    f = lambda x: evaluate(F, T, [x, F(9), F(11)])
    assert f(F(2)) - f(F(1)) == f(F(1)) - f(F(0))


if __name__ == "__main__":
    _path.run_module(globals())
