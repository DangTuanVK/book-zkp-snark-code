"""Kiểm Chương 5: DSL (R1CS / Plonkish / CCS cùng thỏa), Groth16, tái ngẫu nhiên hóa."""
import random
import _path  # noqa: F401
from snarklib.field import GF
from snarklib.curve import toy_generator
from snarklib.r1cs import factor_r1cs
from snarklib.dsl import Circuit, mach_phan_tich
from snarklib.groth16 import (r1cs_to_qap, setup, prove, verify, rerandomize,
                              simulate, combine)

F = GF(449)
E, P, G = toy_generator()
DUNG = {"N": 391, "a": 17, "b": 23}


def _ba_ban(c, inp):
    r, pk, ccs = c.compile_r1cs(), c.compile_plonkish(), c.compile_ccs()
    x, w = c.r1cs_assignment(inp)
    return (r.is_satisfied(x, w), pk.is_satisfied(*c.plonkish_assignment(inp)),
            ccs.is_satisfied(x, w))


def test_ba_ban_bien_dich_cung_thoa():
    for nz in (False, True):
        c = mach_phan_tich(nonzero=nz)[0]
        assert _ba_ban(c, DUNG) == (True, True, True)
        assert _ba_ban(c, {"N": 391, "a": 17, "b": 24}) == (False, False, False)
        assert _ba_ban(c, {"N": 392, "a": 17, "b": 23}) == (False, False, False)
        assert _ba_ban(c, {"N": 391 + 449 * 0, "a": 33, "b": 1}) == (False, False, False)


def test_so_rang_buoc():
    c = mach_phan_tich()[0]
    st = c.stats()
    assert (st["m"], st["n"], st["n_public"]) == (13, 14, 1)
    pk = c.compile_plonkish()
    assert pk.n == 20 and sum(1 for v in pk.sel["qM"] if v) == 11
    assert c.compile_plonkish(pad_pow2=True).n == 32
    assert c.compile_ccs().degree() == 2


def test_tuong_thich_factor_r1cs():
    # Hàng nhân cuối cùng, cắt về 4 cột đầu z = (1, N, a, b), chính là factor_r1cs
    r = mach_phan_tich()[0].compile_r1cs()
    ref = factor_r1cs(F)
    assert [r.A[-1][:4]] == ref.A and [r.B[-1][:4]] == ref.B and [r.C[-1][:4]] == ref.C
    assert all(v == 0 for M in (r.A, r.B, r.C) for v in M[-1][4:])


def test_witness_xau_bi_bac_o_moi_dang():
    # Sửa tay một bit thành 2 (không boolean) nhưng giữ tổng: R1CS và Plonkish đều bác
    c = mach_phan_tich()[0]
    r = c.compile_r1cs()
    x, w = c.r1cs_assignment(DUNG)
    w = list(w)
    # a = 17 = 1 + 16: a0 = 1, a4 = 1 -> thay (a0, a1) = (1, 0) bằng (3, -1)
    w[2], w[3] = F(3), F(-1)
    assert not r.is_satisfied(x, w)
    assert not c.compile_ccs().is_satisfied(x, w)


def test_vet_can_nghiem():
    c = mach_phan_tich()[0]
    r = c.compile_r1cs()
    sols = []
    for a in range(32):
        for b in range(32):
            x, w = c.r1cs_assignment({"N": 391, "a": a, "b": b})
            if r.is_satisfied(x, w):
                sols.append((a, b))
    assert sols == [(17, 23), (23, 17), (28, 30), (30, 28)]


def test_dsl_mul_assert_eq():
    c = Circuit()
    x = c.public("x")
    y = c.private("y")
    y2 = c.mul(y, y)
    c.assert_eq(y2 + y, x)            # y^2 + y = x
    for inp, ok in (({"x": 12, "y": 3}, True), ({"x": 12, "y": 4}, False)):
        assert _ba_ban(c, inp) == (ok, ok, ok)


def _chuan_bi(seed=2016):
    c = mach_phan_tich()[0]
    qap = r1cs_to_qap(c.compile_r1cs())
    pk, vk, tau = setup(qap, G, G, random.Random(seed))
    return c, qap, pk, vk, tau


def test_qap():
    c, qap, _, _, _ = _chuan_bi()
    assert qap.n == 16
    x, w = c.r1cs_assignment(DUNG)
    _, _, _, h, rem = combine(qap, [1] + x + w)
    assert rem.is_zero() and h.degree() <= qap.n - 2
    x, w = c.r1cs_assignment({"N": 391, "a": 17, "b": 24})
    assert not combine(qap, [1] + x + w)[4].is_zero()


def test_groth16_dung_sai():
    c, qap, pk, vk, _ = _chuan_bi()
    x, w = c.r1cs_assignment(DUNG)
    pi = prove(pk, qap, [1] + x + w, random.Random(1))
    assert verify(vk, [391], pi)
    assert not verify(vk, [392], pi)
    assert not verify(vk, [391], pi._replace(C=pi.C + G))
    # Pairing đối xứng (Type 1): đổi chỗ A và B vẫn verify -- ở Type 3 không gõ kiểu được
    assert verify(vk, [391], pi._replace(A=pi.B, B=pi.A))


def test_tai_ngau_nhien():
    c, qap, pk, vk, _ = _chuan_bi()
    x, w = c.r1cs_assignment(DUNG)
    pi = prove(pk, qap, [1] + x + w, random.Random(1))
    for r, s in [(5, 7), (1, 1), (448, 0), (100, 300)]:
        moi = rerandomize(vk, pi, r, s)
        assert verify(vk, [391], moi)
        assert (r, s) == (1, 0) or moi != pi


def test_mo_phong():
    c, qap, pk, vk, tau = _chuan_bi()
    for N in (391, 389):
        assert verify(vk, [N], simulate(qap, tau, [N], G, G, random.Random(N)))


if __name__ == "__main__":
    _path.run_module(globals())
