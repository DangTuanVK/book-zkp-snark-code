"""Kiểm plonkish.py: cổng + ràng buộc sao chép trên ví dụ a·b = 391."""
import _path  # noqa: F401
from snarklib.field import GF
from snarklib.plonkish import factor_plonkish

F = GF(449)


def test_plonkish():
    c = factor_plonkish(F)
    assert c.is_satisfied(a=[391, 17], b=[0, 23], c=[0, 391], x=[391])
    assert not c.is_satisfied(a=[391, 17], b=[0, 23], c=[0, 391], x=[392])
    # cổng đúng nhưng vi phạm sao chép
    assert c.copy_violations(a=[390, 17], b=[0, 23], c=[0, 391])
    assert not c.is_satisfied(a=[391, 17], b=[0, 22], c=[0, 391], x=[391])


if __name__ == "__main__":
    _path.run_module(globals())
