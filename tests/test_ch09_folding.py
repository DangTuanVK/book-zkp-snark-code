"""Test Chương 9: chu trình/2-chain, số học không gốc, folding Nova và HyperNova."""
import random
import _path  # noqa: F401
from snarklib.field import GF
from snarklib.transcript import Transcript
from snarklib.dsl import mach_phan_tich
from snarklib.ccs import CCS
from snarklib import folding as fo

F = GF(449)
rng = random.Random(9)
CIRC = mach_phan_tich()[0]


def _wit(a, b):
    return CIRC.r1cs_assignment({"N": 391, "a": a, "b": b})


def test_count_points_vs_brute():
    for p, a, b in ((7, 0, 3), (13, 0, 6), (31, 2, 5), (43, 1, 1)):
        brute = 1 + sum(1 for x in range(p) for y in range(p) if (y * y - x ** 3 - a * x - b) % p == 0)
        assert fo.count_points(p, a, b) == brute


def test_cycles():
    cyc = fo.find_cycles(5, 400, b_max=10, limit=3)
    assert cyc[0]["q1"] == 7 and cyc[0]["q2"] == 13
    for c in cyc:
        assert fo.count_points(c["q1"], 0, c["b1"]) == c["q2"]
        assert fo.count_points(c["q2"], 0, c["b2"]) == c["q1"]


def test_chain_outer():
    out = fo.chain_outer_supersingular(5387)
    assert out["P"] == 129287 and out["n"] == 24 * 5387 and out["order"] == 5387 and out["k"] == 2
    assert fo.ec_mul(5387, out["G"], 1, out["P"]) is None


def test_nonnative_formula_and_soundness():
    for w in (1, 2):
        c, prm, assign = fo.nonnative_mul_circuit(5387, F, w)
        assert c.stats()["m"] == fo.nonnative_formula(prm["l"], w, prm["wc"])
        for _ in range(5):
            a, b = rng.randrange(5387), rng.randrange(5387)
            assert c.is_satisfied_env(c.solve(assign(a, b)))
            # c sai lệch 1: không thỏa
            assert not c.is_satisfied_env(c.solve(assign(a, b, cc=a * b % 5387 + 1, kk=a * b // 5387)))
    assert not fo.nonnative_params(5387, 449, 3)["ok"]


def test_nova_fold():
    R = fo.RelaxedR1CS(CIRC.compile_r1cs())
    i1, W1 = R.trivial(*_wit(17, 23))
    I, W = i1, W1
    for k, (a, b) in enumerate([(23, 17), (17, 23), (28, 30)]):
        i2, W2 = R.trivial(*_wit(a, b))
        I, W, _, _ = R.fold(I, W, i2, W2, Transcript(b"t%d" % k))
        assert R.is_satisfied(I, W)
    bad = 0
    for k in range(20):
        ib, Wb = R.trivial(*_wit(17, 24))
        bad += R.is_satisfied(*R.fold(i1, W1, ib, Wb, Transcript(b"n%d" % k))[:2])
    assert bad <= 2


def test_hypernova_complete_and_chain():
    P = fo.CCSPadded(CIRC.compile_ccs())
    x1, w1 = _wit(17, 23)
    tp, tv = Transcript(b"h"), Transcript(b"h")
    L, w = P.linearize(x1, w1, tp), w1
    Lv = P.linearize(x1, w1, tv)
    for k, (a, b) in enumerate([(23, 17), (17, 23), (23, 17)]):
        x2, w2 = _wit(a, b)
        pf, L, w = P.hypernova_prove(L, w, x2, w2, tp, label=b"s%d" % k)
        ok, Lv = P.hypernova_verify(Lv, x2, pf, tv, label=b"s%d" % k)
        assert ok and Lv == L and P.lcccs_is_satisfied(L, w)


def test_hypernova_reject():
    P = fo.CCSPadded(CIRC.compile_ccs())
    x1, w1 = _wit(17, 23)
    x2, w2 = _wit(23, 17)
    rej_sigma = rej_v = 0
    for k in range(20):
        lab = b"r%d" % k
        tp, tv = Transcript(lab), Transcript(lab)
        pf, _, _ = P.hypernova_prove(P.linearize(x1, w1, tp), w1, x2, w2, tp)
        pf["sigma"] = [pf["sigma"][0] + 1] + pf["sigma"][1:]
        rej_sigma += not P.hypernova_verify(P.linearize(x1, w1, tv), x2, pf, tv)[0]
        # LCCCS khai man (v_1 + 1): Verifier hoặc từ chối, hoặc LCCCS mới không thỏa
        tp = Transcript(lab)
        L1 = P.linearize(x1, w1, tp)
        L1["v"] = [L1["v"][0] + 1] + L1["v"][1:]
        pf, L, w = P.hypernova_prove(L1, w1, x2, w2, tp)
        tv = Transcript(lab)
        Lb = P.linearize(x1, w1, tv)
        Lb["v"] = list(L1["v"])
        ok, _ = P.hypernova_verify(Lb, x2, pf, tv)
        rej_v += not (ok and P.lcccs_is_satisfied(L, w))
    assert rej_sigma >= 18 and rej_v >= 18


def test_hypernova_degree3():
    # CCS bậc 3: z = (1, d, e, a, b, c); a·b·c = d, b·c·c = e
    def sel(i):
        return [[1 if j == i else 0 for j in range(6)]]
    M = [sel(3) + sel(4), sel(4) + sel(5), sel(5) + sel(5), sel(1) + sel(2)]
    ccs = CCS(F, M, [[0, 1, 2], [3]], [1, -1], n_public=2)
    P = fo.CCSPadded(ccs)
    assert P.d == 3 and P.s == 1
    a, b, cc = 5, 7, 11
    x1, w1 = [a * b * cc % 449, b * cc * cc % 449], [a, b, cc]
    x2, w2 = [2 * 3 * 4, 3 * 16], [2, 3, 4]
    assert ccs.is_satisfied(x1, w1) and ccs.is_satisfied(x2, w2)
    tp, tv = Transcript(b"d3"), Transcript(b"d3")
    pf, L, w = P.hypernova_prove(P.linearize(x1, w1, tp), w1, x2, w2, tp)
    assert len(pf["sc"]["polys"][0]) == 5          # bậc d + 1 = 4 -> 5 giá trị
    ok, Lv = P.hypernova_verify(P.linearize(x1, w1, tv), x2, pf, tv)
    assert ok and Lv == L and P.lcccs_is_satisfied(L, w)


if __name__ == "__main__":
    _path.run_module(globals())
