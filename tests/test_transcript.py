"""Kiểm transcript.py: tất định, nhạy với nhãn/thứ tự/dữ liệu."""
import _path  # noqa: F401
from snarklib.transcript import Transcript
from snarklib.field import GF
from snarklib.curve import toy_generator


def run(label, items):
    t = Transcript(label)
    for lab, d in items:
        t.absorb(lab, d)
    return t.challenge("alpha"), t.challenge("beta")


def test_deterministic():
    E, P, G = toy_generator()
    F = GF(449)
    items = [("N", 391), ("C", G), ("v", F(17)), ("raw", b"\x01\x02")]
    a1 = run("vi-du", items)
    a2 = run("vi-du", items)
    assert a1 == a2
    assert a1[0] != a1[1]                     # hai thách thức liên tiếp khác nhau
    assert isinstance(a1[0], F) and 0 <= int(a1[0]) < 449


def test_sensitivity():
    base = run("vi-du", [("N", 391), ("x", 5)])
    assert run("vi-du-2", [("N", 391), ("x", 5)]) != base       # nhãn giao thức
    assert run("vi-du", [("x", 5), ("N", 391)]) != base         # thứ tự
    assert run("vi-du", [("N", 391), ("x", 6)]) != base         # dữ liệu
    assert run("vi-du", [("M", 391), ("x", 5)]) != base         # nhãn absorb
    assert run("vi-du", [("N", b"\x01\x87"), ("x", 5)]) != base  # thẻ kiểu


def test_distribution():
    # 4490 thách thức: tần suất mỗi giá trị trung bình 10, không giá trị nào lệch quá xa
    F = GF(449)
    counts = [0] * 449
    t = Transcript("phan-bo")
    for k in range(4490):
        t.absorb("k", k)
        counts[int(t.challenge("c"))] += 1
    assert max(counts) < 40 and sum(counts) == 4490


if __name__ == "__main__":
    _path.run_module(globals())
