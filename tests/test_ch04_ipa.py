"""Kiểm pedersen.py và ipa.py (Chương 4): đúng, sai, proof hỏng bị bác."""
import _path  # noqa: F401
from snarklib.field import GF
from snarklib.transcript import Transcript
from snarklib.pedersen import (PedersenParams, hash_to_curve, commit, dlog_bruteforce,
                               E, msm, padd)
from snarklib import ipa
from snarklib.curve import toy_generator

F = GF(449)


def test_pedersen_generators():
    pp = PedersenParams(8, b"test")
    pts = pp.G + [pp.H, pp.U]
    assert len(set(pts)) == len(pts)                       # đôi một khác nhau
    for P in pts:
        assert not P.is_infinity() and (449 * P).is_infinity()
    assert hash_to_curve(b"test/G", 3) == hash_to_curve(b"test/G", 3)   # tất định
    # đồng cấu: Com(a) + Com(b) = Com(a + b)
    a, b = [1, 2, 3], [5, 7, 11]
    assert padd(commit(pp.G, a), commit(pp.G, b)) == commit(pp.G, [x + y for x, y in zip(a, b)])
    # trên toy, log rời rạc giải được -> phá binding (chỉ minh họa)
    _, _, G1 = toy_generator()
    k = dlog_bruteforce(G1, pp.G[0])
    assert k * G1 == pp.G[0]


def test_ipa_pcs():
    for n in (2, 8, 32):
        pp = ipa.setup(n)
        a = [(3 * i * i + 1) % 449 for i in range(n)]
        P = ipa.commit(pp, a)
        z = 17
        v, pr = ipa.open(pp, a, z, Transcript(b"t"))
        assert v == sum(F(c) * F(z) ** i for i, c in enumerate(a))
        assert len(pr["L"]) == n.bit_length() - 1
        assert ipa.verify(pp, P, z, v, pr, Transcript(b"t"))
        assert not ipa.verify(pp, P, z, v + 1, pr, Transcript(b"t"))      # sai giá trị
        assert not ipa.verify(pp, P, z + 1, v, pr, Transcript(b"t"))      # sai điểm
        assert not ipa.verify(pp, P, z, v, pr, Transcript(b"khac"))      # transcript khác
        bad = dict(pr, a=pr["a"] + 1)
        assert not ipa.verify(pp, P, z, v, bad, Transcript(b"t"))
        if n > 1:
            bad = dict(pr, L=[pr["L"][0] + pp.G[0]] + pr["L"][1:])
            assert not ipa.verify(pp, P, z, v, bad, Transcript(b"t"))


def test_b_final_and_s():
    us = [F(5), F(7), F(11)]
    s = ipa.s_vector(us)
    assert len(s) == 8 and s[0] == (F(5) * 7 * 11).inverse() and s[7] == F(5 * 7 * 11)
    z = 13
    assert ipa.b_final_powers(us, z) == ipa.inner(s, ipa.powers(z, 8))


def test_bulletproofs_general():
    n = 8
    G = PedersenParams(n, b"bpG").G
    H = PedersenParams(n, b"bpH").G
    U = hash_to_curve(b"bpU", 0)[0]
    a = [2, 3, 5, 7, 11, 13, 17, 19]
    b = [1, 4, 9, 16, 25, 36, 49, 64]
    c = ipa.inner(a, b)
    P = padd(msm(a, G), msm(b, H))
    pr = ipa.bp_prove(G, H, U, a, b, Transcript(b"bp"))
    for m in ("fold", "msm"):
        assert ipa.bp_verify(G, H, U, P, c, pr, Transcript(b"bp"), m)
        assert not ipa.bp_verify(G, H, U, P, c + 1, pr, Transcript(b"bp"), m)
        assert not ipa.bp_verify(G, H, U, P, c, dict(pr, b=pr["b"] + 1), Transcript(b"bp"), m)


def test_one_round_identity():
    """Sau một vòng: x^2 L + P + x^-2 R = <a',G'> + <a',b'>U (đẳng thức gập đôi)."""
    pp = ipa.setup(4)
    a = [F(v) for v in (3, 1, 4, 1)]
    b = ipa.powers(5, 4)
    P = padd(msm(a, pp.G), ipa.inner(a, b) * pp.U)
    L, R, _, _ = ipa.round_messages(pp.G, a, b, pp.U)
    x = F(10)
    G2, a2, b2 = ipa.fold_vectors(pp.G, a, b, x)
    assert ipa.fold_commitment(P, L, R, x) == padd(msm(a2, G2), ipa.inner(a2, b2) * pp.U)
    assert E.infinity() == 0 * pp.U


if __name__ == "__main__":
    _path.run_module(globals())
