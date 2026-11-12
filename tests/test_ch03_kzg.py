"""Kiểm Chương 3: KZG (một điểm, batch, đa điểm, GWC, SHPLONK, hiding), PST, ceremony, opcount."""
import importlib.util
import os
import random

import _path  # noqa: F401
from snarklib import kzg, pst, opcount as oc
from snarklib.poly import Poly
from snarklib.mle import evaluate

F = kzg.F
rng = random.Random(2026)
SRS = kzg.setup(15, 77, g2_degree=4)


def _rand_poly(d=15):
    return Poly(F, [rng.randrange(449) for _ in range(d + 1)])


def _tr(label=b"test-ch03"):
    return kzg.fresh_transcript(label)


def test_single_point():
    for _ in range(5):
        f = _rand_poly()
        C = kzg.commit(SRS, f)
        z = rng.randrange(449)
        y, W = kzg.open(SRS, f, z)
        assert y == f(F(z)) and kzg.verify(SRS, C, z, y, W)
        assert not kzg.verify(SRS, C, z, y + 1, W)
        assert not kzg.verify(SRS, C, (z + 1) % 449, y, W)
        assert not kzg.verify(SRS, kzg.commit(SRS, f + 1), z, y, W)


def test_commit_homomorphic_and_degree_bound():
    f, g = _rand_poly(), _rand_poly()
    assert kzg.commit(SRS, f + g) == kzg.commit(SRS, f) + kzg.commit(SRS, g)
    try:
        kzg.commit(SRS, _rand_poly(16))
        assert False, "phải báo lỗi vượt bậc SRS"
    except ValueError:
        pass


def test_batch_one_point():
    fs = [_rand_poly() for _ in range(5)]
    Cs = [kzg.commit(SRS, f) for f in fs]
    ys, W = kzg.open_batch(SRS, fs, Cs, 11, _tr())
    assert kzg.verify_batch(SRS, Cs, 11, ys, W, _tr())
    for i in range(5):
        bad = list(ys)
        bad[i] = bad[i] + 1
        assert not kzg.verify_batch(SRS, Cs, 11, bad, W, _tr())
    assert not kzg.verify_batch(SRS, Cs, 11, ys, W, _tr(b"nhan-khac"))


def test_multi_point_one_poly():
    f = _rand_poly()
    C = kzg.commit(SRS, f)
    for S in ([1, 2], [3, 5, 7], [0, 100, 200, 300]):
        ys, W = kzg.open_multi(SRS, f, S)
        assert kzg.verify_multi(SRS, C, S, ys, W)
        bad = list(ys)
        bad[-1] = bad[-1] + 1
        assert not kzg.verify_multi(SRS, C, S, bad, W)


def test_hand_example():
    f = Poly(F, [23, 17, 0, 1])
    ys = [f(F(1)), f(F(2))]
    assert [int(y) for y in ys] == [41, 65]
    psi, rem = (f - Poly(F, [17, 24])).divmod(kzg.vanishing([1, 2]))
    assert rem.is_zero() and psi == Poly(F, [3, 1])


SETS = [[[3, 7], [3], [7, 11], [3, 7]], [[5]] * 3, [[1, 2, 3], [4], [2, 4]]]


def _rejected_mostly(results):
    """Trường toy r = 449: một bằng chứng sai lọt qua với xác suất ~1/449 cho mỗi
    thách thức, nên phép kiểm âm chạy với 4 nhãn transcript và cho phép lọt tối đa 1."""
    return sum(1 for ok in results if ok) <= 1


LABELS = [b"neg-0", b"neg-1", b"neg-2", b"neg-3"]


def test_gwc_multi_poly_multi_point():
    for Ss in SETS:
        fs = [_rand_poly() for _ in Ss]
        Cs = [kzg.commit(SRS, f) for f in fs]
        pr = kzg.open_gwc(SRS, fs, Cs, Ss, _tr())
        assert len(pr.Ws) == len({p for S in Ss for p in S})
        assert kzg.verify_gwc(SRS, Cs, Ss, pr, _tr())
        res = []
        for lab in LABELS:
            pr = kzg.open_gwc(SRS, fs, Cs, Ss, _tr(lab))
            bad = [list(y) for y in pr.ys]
            bad[0][0] = bad[0][0] + 1
            res.append(kzg.verify_gwc(SRS, Cs, Ss, pr._replace(ys=bad), _tr(lab)))
        assert _rejected_mostly(res)


def test_shplonk():
    for Ss in SETS:
        fs = [_rand_poly() for _ in Ss]
        Cs = [kzg.commit(SRS, f) for f in fs]
        pr = kzg.open_shplonk(SRS, fs, Cs, Ss, _tr())
        with oc.scope() as c:
            assert kzg.verify_shplonk(SRS, Cs, Ss, pr, _tr())
        assert c["pairing"] == 2
        Cs2 = list(Cs)
        Cs2[0] = Cs2[0] + SRS.G1
        for i in range(len(Ss)):
            res = []
            for lab in LABELS:
                pr = kzg.open_shplonk(SRS, fs, Cs, Ss, _tr(lab))
                bad = [list(y) for y in pr.ys]
                bad[i][-1] = bad[i][-1] + 1
                res.append(kzg.verify_shplonk(SRS, Cs, Ss, pr._replace(ys=bad), _tr(lab)))
            assert _rejected_mostly(res)
        res = []
        for lab in LABELS:
            pr = kzg.open_shplonk(SRS, fs, Cs, Ss, _tr(lab))
            # lưu ý: khi T chỉ có một điểm, W' = W (L(X) = (X - z)h(X)) nên không thử hoán đổi
            res.append(kzg.verify_shplonk(SRS, Cs, Ss, pr._replace(W=pr.W + SRS.G1), _tr(lab)))
            res.append(kzg.verify_shplonk(SRS, Cs2, Ss, pr, _tr(lab)))
        assert _rejected_mostly(res)


def test_hiding():
    hs = kzg.hiding_setup(15, 77)
    f, fh = _rand_poly(), _rand_poly()
    C = kzg.commit_hiding(hs, f, fh)
    assert C != kzg.commit(hs, f)
    y, yh, W = kzg.open_hiding(hs, f, fh, 9)
    assert kzg.verify_hiding(hs, C, 9, y, yh, W)
    assert not kzg.verify_hiding(hs, C, 9, y + 1, yh, W)
    H, _ = kzg.hash_to_g1()
    beta = kzg.dlog_bruteforce(hs.G1, H)
    assert beta is not None and beta * hs.G1 == H        # hạn chế của toy


def test_pst():
    for l in (1, 2, 3, 4):
        taus = [rng.randrange(1, 449) for _ in range(l)]
        srs = pst.setup(l, taus)
        T = [rng.randrange(449) for _ in range(2 ** l)]
        C = pst.commit(srs, T)
        z = [rng.randrange(449) for _ in range(l)]
        v, Ws = pst.open(srs, T, z)
        assert v == evaluate(F, T, z) and len(Ws) == l
        assert pst.verify(srs, C, z, v, Ws)
        assert not pst.verify(srs, C, z, v + 1, Ws)
        # commit = [q~(τ)]_1
        assert C == int(evaluate(F, T, taus)) * srs.G1


def _ceremony_module():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(root, "ch03", "ma_nguon_03_03_ceremony.py")
    spec = importlib.util.spec_from_file_location("ceremony_ch03", path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_ceremony_detects_cheating():
    m = _ceremony_module()
    srs0 = {"g1": [m.G1] * (m.D + 1), "g2": [m.G2, m.G2]}
    s1, bc1 = m.dong_gop(srs0, 17)
    assert m.kiem_dong_gop(srs0, s1, bc1) == []
    s2, bc2 = m.dong_gop(s1, 23)
    assert m.kiem_dong_gop(s1, s2, bc2) == []
    assert s2["g1"][3] == pow(17 * 23, 3, 449) * m.G1
    # thay một phần tử -> K5 tại i = 4 và 5
    x, b = m.thay_mot_phan_tu(s1, s2, bc2)
    loi = m.kiem_dong_gop(s1, x, b)
    assert loi and all(e.startswith("K5") for e in loi) and "i = 4" in loi[0]
    # bỏ lịch sử -> K2
    x, b = m.bo_lich_su(s1, s2, bc2)
    assert any(e.startswith("K2") for e in m.kiem_dong_gop(s1, x, b))
    # [s]_1, [s]_2 lệch nhau -> K3
    b = dict(bc2, s1=bc2["s1"] + m.G1)
    assert any(e.startswith("K3") for e in m.kiem_dong_gop(s1, s2, b))
    # [τ]_2 không khớp [τ]_1 -> K4 (và K5)
    x = {"g1": s2["g1"], "g2": [m.G2, s2["g2"][1] + m.G2]}
    assert any(e.startswith("K4") for e in m.kiem_dong_gop(s1, x, bc2))


def test_opcount_scope():
    G = SRS.G1
    with oc.scope() as outer:
        with oc.scope() as inner:
            oc.msm([G, G, G], [0, 1, 5])
            oc.smul(G, 6)
        oc.pairing_check([(G, G), (-G, G)])
    assert inner["msm"] == 1 and inner["msm_size"] == 3 and inner["msm_nontrivial"] == 1
    assert inner["smul"] == 1 and "pairing" not in inner
    assert outer["msm"] == 1 and outer["pairing"] == 2 and outer["final_exp"] == 1
    assert oc.smul(G, 448) == -G


if __name__ == "__main__":
    _path.run_module(globals())
