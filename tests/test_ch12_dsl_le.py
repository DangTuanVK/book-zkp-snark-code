"""Kiểm gadget assert_le và mạch đã sửa (thêm vào snarklib.dsl ở Chương 12)."""
import _path  # noqa: F401
from snarklib.dsl import Circuit, mach_phan_tich, mach_phan_tich_sua


def _thoa(c, inp):
    x, w = c.r1cs_assignment(inp)
    return c.compile_r1cs().is_satisfied(x, w)


def test_assert_le_dung_mien():
    for B in (0, 1, 7, 19, 23, 31):
        c = Circuit()
        x = c.private("x")
        c.to_bits(x, 5, prefix="x")
        bits = c.assert_le(x, B)
        assert len(bits) == max(1, B.bit_length())
        ok = [v for v in range(32) if _thoa(c, {"x": v})]
        assert ok == list(range(B + 1)), (B, ok)


def test_assert_le_so_rang_buoc():
    c = Circuit()
    x = c.private("x")
    c.assert_le(x, 19)
    st = c.stats()
    assert st["m"] == 6 and st["labels"] == {"bool": 5, "le": 1}


def test_assert_le_tham_so_sai():
    c = Circuit()
    x = c.private("x")
    for bad in ((-1, None), (19, 4)):
        try:
            c.assert_le(x, bad[0], k=bad[1])
            assert False
        except ValueError:
            pass


def test_mach_cu_khong_doi():
    c = mach_phan_tich()[0]
    assert c.stats()["m"] == 13 and c.compile_plonkish().n == 20


def test_mach_da_sua():
    c = mach_phan_tich_sua()[0]
    assert c.stats()["m"] == 25 and c.compile_plonkish().n == 38
    assert _thoa(c, {"N": 391, "a": 17, "b": 23})
    for a, b in ((23, 17), (28, 30), (30, 28), (17, 24)):
        assert not _thoa(c, {"N": 391, "a": a, "b": b})
    # ba bản biên dịch cùng thỏa
    inp = {"N": 391, "a": 17, "b": 23}
    x, w = c.r1cs_assignment(inp)
    assert c.compile_plonkish().is_satisfied(*c.plonkish_assignment(inp))
    assert c.compile_ccs().is_satisfied(x, w)


if __name__ == "__main__":
    _path.run_module(globals())
