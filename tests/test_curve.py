"""Kiểm curve.py: luật nhóm, đếm điểm, sinh G1 tất định."""
import _path  # noqa: F401
from snarklib.curve import toy_curve, toy_generator


def test_count():
    E = toy_curve()
    assert E.count_points() == 5388 == 12 * 449


def test_group_law():
    E, P, G = toy_generator()
    O = E.infinity()
    assert P.on_curve() and G.on_curve()
    assert G + O == G and G - G == O
    A, B = 5 * G, 11 * G
    assert A + B == B + A == 16 * G
    assert (A + B) + G == A + (B + G)
    assert 2 * G == G + G and (-3) * G == -(3 * G)
    assert (5388 * P).is_infinity()


def test_generator():
    E, P, G = toy_generator()
    # x nhỏ nhất thỏa điều kiện; mọi x nhỏ hơn đều bị loại
    for x in range(int(P.x)):
        v = E.rhs(x)
        assert v == 0 or v.legendre() != 1 or (12 * E.point(x, v.sqrt())).is_infinity()
    assert int(P.y) < 5387 - int(P.y)
    assert not G.is_infinity() and (449 * G).is_infinity()
    assert G == 12 * P and G.order(449) == 449


if __name__ == "__main__":
    _path.run_module(globals())
