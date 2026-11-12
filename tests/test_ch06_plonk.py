"""Test mini-PLONK (Chương 6): snarklib/plonk.py."""
import random

import _path  # noqa: F401
from snarklib import plonk, kzg
from snarklib.plonk import F
from snarklib.dsl import mach_phan_tich

_c = mach_phan_tich()[0]
CIRC = _c.compile_plonkish()
A, B, C, X = _c.plonkish_assignment({"N": 391, "a": 17, "b": 23})
PK, VK = plonk.setup(CIRC, tau=123)
LABELS = [b"t6/1", b"t6/2", b"t6/3", b"t6/4"]


def _honest(seed=1, label=b"snarklib/plonk"):
    return plonk.prove(PK, A, B, C, X, random.Random(seed), label=label)[0]


def test_tham_so_mien():
    assert PK.n == 32 and int(PK.omega) == 221 and int(PK.k1) == 2 and int(PK.k2) == 3
    assert len(PK.srs.g1) == PK.n + 6          # [τ^0..τ^{n+5}]_1


def test_proof_trung_thuc():
    for s in range(3):
        assert plonk.verify(VK, X, _honest(s))


def test_coset_fft_bang_chia_da_thuc():
    p = plonk.Prover(PK, A, B, C, X, random.Random(4))
    p.run()
    t, rem = p.quotient_naive()
    assert rem.is_zero() and p.quotient_coset_fft() == t and t.degree() <= 3 * PK.n + 5


def test_z_va_mo_kzg_rieng():
    p = plonk.Prover(PK, A, B, C, X, random.Random(5))
    pf = p.run()
    assert p.z_vals[0] == 1
    # W_ζω là nhân chứng KZG một điểm của z tại ζω
    assert kzg.verify(PK.srs, pf.z, p.zeta * PK.omega, pf.z_omega_bar, pf.W_zeta_omega)


def test_sua_day_bi_bac():
    a2 = list(A)
    a2[19] = F(18)
    lot = 0
    for lab in LABELS:
        pf, _ = plonk.prove(PK, a2, B, C, X, random.Random(7), label=lab, check=False)
        lot += plonk.verify(VK, X, pf, label=lab)
    assert lot <= 1
    try:
        plonk.prove(PK, a2, B, C, X, random.Random(7))
        assert False, "Prover trung thực phải từ chối witness sai"
    except ValueError:
        pass


def test_doi_dau_vao_cong_khai_bi_bac():
    lot = sum(plonk.verify(VK, [389], _honest(8, lab), label=lab) for lab in LABELS)
    assert lot <= 1


def test_sua_tung_phan_tu_proof():
    pf = _honest(9)
    for k in plonk.FR_FIELDS:
        assert not plonk.verify(VK, X, pf._replace(**{k: getattr(pf, k) + 1}))
    for k in plonk.G1_FIELDS:
        bad = pf._replace(**{k: getattr(pf, k) + PK.srs.G1})
        assert not plonk.verify(VK, X, bad)
    assert not plonk.verify(VK, [391, 1], pf)


def test_setup_tu_choi_tau_trong_H():
    try:
        plonk.setup(CIRC, tau=77)               # 77^32 = 1 trong F_449
        assert False
    except ValueError:
        pass


def test_frozen_heart_mot_lan_thu():
    pf, xs, lan = plonk.forge_frozen_heart(PK, A, B, C, [391], None, random.Random(3))
    assert lan == 1 and xs != 391
    assert plonk.verify(VK, [xs], pf, hash_public=False)
    assert not plonk.verify(VK, [xs], pf, hash_public=True)


def test_frozen_heart_nham_389():
    pf, xs, lan = plonk.forge_frozen_heart(PK, A, B, C, [391], 389, random.Random(3))
    assert xs == 389 and plonk.verify(VK, [389], pf, hash_public=False)
    assert not plonk.verify(VK, [389], pf, hash_public=True)


if __name__ == "__main__":
    _path.run_module(globals())
