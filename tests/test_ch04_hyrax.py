"""Kiểm hyrax.py: mở MLE 16 hệ số (ma trận 4 × 4), cả hai chế độ, proof hỏng bị bác."""
import _path  # noqa: F401
from snarklib.field import GF
from snarklib.transcript import Transcript
from snarklib import hyrax, mle

F = GF(449)


def test_hyrax_16():
    T = [(5 * i + 2) % 449 for i in range(16)]
    pp = hyrax.setup(16)
    comms, st = hyrax.commit(pp, T)
    assert len(comms) == 4 and pp.n == 4
    r = [F(2), F(3), F(5), F(7)]
    for mode in ("ipa", "truc_tiep"):
        v, pr = hyrax.open(pp, st, r, Transcript(b"h"), mode)
        assert v == mle.evaluate(F, T, r)
        assert hyrax.verify(pp, comms, r, v, pr, Transcript(b"h"), 16)
        assert not hyrax.verify(pp, comms, r, v + 1, pr, Transcript(b"h"), 16)
        bad = list(comms)
        bad[1] = bad[1] + pp.G[0]                     # cam kết hàng bị sửa
        assert not hyrax.verify(pp, bad, r, v, pr, Transcript(b"h"), 16)
    assert hyrax.sizes(16) == {"commit_group": 4, "proof_group": 4, "proof_field": 1}


def test_hyrax_odd_vars():
    T = list(range(1, 33))                             # l = 5: ma trận 4 × 8
    pp = hyrax.setup(32)
    comms, st = hyrax.commit(pp, T)
    r = [F(v) for v in (9, 8, 7, 6, 5)]
    v, pr = hyrax.open(pp, st, r, Transcript(b"h"))
    assert len(comms) == 4 and v == mle.evaluate(F, T, r)
    assert hyrax.verify(pp, comms, r, v, pr, Transcript(b"h"), 32)


if __name__ == "__main__":
    _path.run_module(globals())
