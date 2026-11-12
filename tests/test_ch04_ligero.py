"""Kiểm ligero.py: mở đúng, giá trị sai, cột/đường Merkle/u_prox hỏng đều bị bác."""
import _path  # noqa: F401
from snarklib.field import GF
from snarklib.transcript import Transcript
from snarklib import ligero, mle

F = GF(449)


def _setup(n):
    T = [(11 * i * i + 3) % 449 for i in range(n)]
    l = n.bit_length() - 1
    r = [F(4 + 3 * i) for i in range(l)]
    root, st = ligero.commit(T)
    pr = ligero.open(st, r, Transcript(b"lg"))
    return T, r, root, st, pr


def test_ligero_ok():
    for n in (4, 8, 16, 64, 256):
        T, r, root, st, pr = _setup(n)
        assert pr["v"] == mle.evaluate(F, T, r)
        assert ligero.verify(root, r, pr["v"], pr, Transcript(b"lg"))
        assert not ligero.verify(root, r, pr["v"] + 1, pr, Transcript(b"lg"))


def test_rs_linear():
    a, b = [1, 2, 3, 4], [7, 0, 5, 448]
    ea, eb = ligero.rs_encode(a), ligero.rs_encode(b)
    eab = ligero.rs_encode([x + 3 * y for x, y in zip(a, b)])
    assert eab == [x + 3 * y for x, y in zip(ea, eb)] and len(ea) == 16


def test_ligero_tampered():
    T, r, root, st, pr = _setup(16)
    v = pr["v"]
    j, col, path = pr["cols"][0]
    bad = dict(pr, cols=[(j, [col[0] + 1] + col[1:], path)] + pr["cols"][1:])
    assert not ligero.verify(root, r, v, bad, Transcript(b"lg"))        # cột bị sửa
    bad = dict(pr, cols=[(j, col, [b"\x00" * 32] + path[1:])] + pr["cols"][1:])
    assert not ligero.verify(root, r, v, bad, Transcript(b"lg"))        # đường Merkle sai
    bad = dict(pr, u_prox=[pr["u_prox"][0] + 1] + pr["u_prox"][1:])
    assert not ligero.verify(root, r, v, bad, Transcript(b"lg"))        # u_prox sai
    # Prover gian: đổi u_eval để khớp giá trị sai -> kiểm cột bắt được
    R = mle.eq_table(F, r[2:])
    ue = list(pr["u_eval"])
    ue[0] = ue[0] + R[0].inverse()                      # <ue, R> tăng đúng 1
    bad = dict(pr, u_eval=ue)
    assert not ligero.verify(root, r, v + 1, bad, Transcript(b"lg"))
    # gốc khác
    root2, _ = ligero.commit([x + 1 for x in T])
    assert not ligero.verify(root2, r, v, pr, Transcript(b"lg"))


def test_sizes():
    s = ligero.sizes(256, t=8)
    assert (s["m"], s["k"], s["ncols"]) == (16, 16, 64)
    assert s["proof_field"] == 2 * 16 + 8 * 16 and s["proof_hash"] == 8 * 6


if __name__ == "__main__":
    _path.run_module(globals())
