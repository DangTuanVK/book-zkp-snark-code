"""Test Chương 7: snarklib/multiset.py, snarklib/logup.py."""
import random

import _path  # noqa: F401
from snarklib import multiset as ms, logup, opcount as oc
from snarklib.multiset import F449
from snarklib.logup import F

LABELS = [b"t7/%d" % i for i in range(12)]


def test_da_tap_tich():
    t = list(range(32))
    m, out = ms.multiplicities([17, 23, 17], t)
    assert out == [] and m[17] == 2 and m[23] == 1 and sum(m) == 3
    assert ms.multiset_equal_poly(F449, [17, 23, 17], ms.expand(t, m))
    assert ms.multiset_check(F449, [1, 2, 3], [3, 1, 2], 77)
    assert not ms.multiset_equal_poly(F449, [17, 40], [17, 23])
    hits, tot = ms.false_positive_rate(F449, [17, 40], [17, 23])
    assert tot == 449 and hits <= 2
    Z = ms.grand_product(F449, [17, 23], [23, 17], 5)
    assert Z[-1] == 1


def test_dong_nhat_thuc_logup():
    t = logup.range_table(31)
    f = [F(17), F(23), F(17)]
    m, out = logup.field_multiplicities(f, t)
    assert out == [] and int(m[17]) == 2
    assert logup.pole_coefficients(f, t, m) == {}
    for lam in (5, 100, 300):
        assert logup.identity_at(f, t, m, lam)
    f2 = [F(17), F(40)]
    m2, out2 = logup.field_multiplicities(f2, t)
    assert out2 == [40] and logup.pole_coefficients(f2, t, m2) == {40: F(1)}
    assert not logup.identity_at(f2, t, m2, 5)


def test_tong_tich_luy_quay_ve_0():
    n, t = 32, logup.range_table(31)
    q = [1 if i in (1, 2) else 0 for i in range(n)]
    f = [F(391), F(17), F(23)] + [F(0)] * (n - 3)
    m, _ = logup.field_multiplicities(f, t, q)
    phi = logup.running_sum(f, t, m, 7, q)
    assert phi[0] == 0 and phi[n] == 0
    f[2] = F(40)
    m, _ = logup.field_multiplicities(f, t, q)
    assert logup.running_sum(f, t, m, 7, q)[n] != 0


def test_dieu_kien_dac_so():
    kq = logup.characteristic_counterexample()
    assert kq["f_len"] == 451 and kq["logup_ok"] and not kq["product_ok_any"]
    assert logup.min_cheat_size() == 449


def test_range_check_can_tuy_y():
    assert logup.range_check([17, 19, 0], 19)[0]
    assert not logup.range_check([17, 20], 19)[0]
    assert not logup.range_check([448], 19)[0]
    assert logup.range_check([23], 23)[0] and not logup.range_check([24], 23)[0]


def test_vet_can():
    assert logup.vet_can(391, (31, 31)) == [(17, 23), (23, 17), (28, 30), (30, 28)]
    assert logup.vet_can(391, (19, 23)) == [(17, 23)]
    assert logup.accepted(391, 28, 30) and not logup.accepted(391, 28, 30, (19, 23))


def test_plonk_logup_trung_thuc():
    for bounds in [(31, 31), (19, 23)]:
        lc = logup.mach_phan_tich_lookup(bounds)
        pk, vk = logup.setup(lc, 123)
        A, B, C, x = logup.assignment(391, 17, 23)
        for s in range(3):
            pf, _ = logup.prove(pk, A, B, C, x, random.Random(s))
            assert logup.verify(vk, x, pf)
        assert logup.proof_size(pf) == (9 + 2 * len(lc.lookups), 6 + 5 * len(lc.lookups))
    assert logup.prove_factor(391, 28, 30, (31, 31))[0]      # nghiệm giả vẫn lọt bảng 5 bit


def test_coset_fft_bang_chia_da_thuc():
    lc = logup.mach_phan_tich_lookup((19, 23))
    pk, _ = logup.setup(lc, 123)
    p = logup.LogupProver(pk, *logup.assignment(391, 17, 23), rng=random.Random(4))
    p.run()
    t, rem = p.quotient_naive()
    assert rem.is_zero() and p.quotient_coset_fft() == t and t.degree() <= 3 * p.n + 5


def test_ngoai_bang_bi_bac():
    lc = logup.mach_phan_tich_lookup((19, 23))
    pk, vk = logup.setup(lc, 123)
    A, B, C, x = logup.assignment(391, 28, 30)
    try:
        logup.prove(pk, A, B, C, x, random.Random(1))
        assert False, "Prover trung thực phải dừng"
    except ValueError:
        pass
    lot = 0
    for i, lab in enumerate(LABELS):
        pf, _ = logup.prove(pk, A, B, C, x, random.Random(50 + i), label=lab, check=False)
        lot += logup.verify(vk, x, pf, label=lab)
    assert lot <= 3          # trên F_449 lọt chủ yếu khi ζ ∈ H (≈ 7%)


def test_sua_proof_bi_bac():
    lc = logup.mach_phan_tich_lookup((31, 31))
    pk, vk = logup.setup(lc, 123)
    A, B, C, x = logup.assignment(391, 17, 23)
    lot = 0
    for i, lab in enumerate(LABELS[:6]):
        pf, _ = logup.prove(pk, A, B, C, x, random.Random(i), label=lab)
        L = pf.lk[0]
        for bad in (L._replace(m_bar=L.m_bar + 1), L._replace(phi_omega_bar=L.phi_omega_bar + 1),
                    L._replace(t_bar=L.t_bar + 1)):
            lot += logup.verify(vk, x, pf._replace(lk=[bad]), label=lab)
        lot += logup.verify(vk, [389], pf, label=lab)
    assert lot <= 1


def test_dem_chi_phi_verifier():
    lc = logup.mach_phan_tich_lookup((19, 23))
    pk, vk = logup.setup(lc, 123)
    A, B, C, x = logup.assignment(391, 17, 23)
    pf, _ = logup.prove(pk, A, B, C, x, random.Random(4))
    with oc.scope() as c:
        assert logup.verify(vk, x, pf)
    assert c["pairing"] == 2 and c["msm_sizes"] == [26, 2]


if __name__ == "__main__":
    _path.run_module(globals())
