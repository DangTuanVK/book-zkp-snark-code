"""Kiểm snarklib.uniqueness: vét cạn có lan truyền, z3 (nếu có), fuzz."""
import random
import _path  # noqa: F401
from snarklib.field import GF
from snarklib.r1cs import factor_r1cs
from snarklib.dsl import mach_phan_tich, mach_phan_tich_sua
from snarklib import uniqueness as U

F = GF(449)


def test_vet_can_mach_cu_va_moi():
    c = mach_phan_tich()[0]
    st = {}
    ws = U.enumerate_witnesses(c.compile_r1cs(), [391], stats=st)
    assert [tuple(w[:2]) for w in ws] == [(17, 23), (23, 17), (28, 30), (30, 28)]
    assert st["complete"] and st["nodes"] == 993
    ws = U.enumerate_witnesses(mach_phan_tich_sua()[0].compile_r1cs(), [391])
    assert [tuple(w[:2]) for w in ws] == [(17, 23)]


def test_vet_can_khop_brute_force():
    c = mach_phan_tich()[0]
    r = c.compile_r1cs()
    for N in (391, 20, 450 - 449, 0):
        ws = sorted(tuple(w[:2]) for w in U.enumerate_witnesses(r, [N]))
        bf = sorted((d["a"], d["b"]) for d in
                    U.brute_force_inputs(c, {"N": N}, {"a": range(32), "b": range(32)}, r))
        assert ws == bf, N


def test_vet_can_ca_truong():
    # a·b = N không kiểm khoảng: 448 nghiệm (a ≠ 0), rẽ nhánh trên cả F_449
    ws = U.enumerate_witnesses(factor_r1cs(F), [391])
    assert len(ws) == 448 and all(a * b % 449 == 391 for a, b in ws)
    assert U.enumerate_witnesses(factor_r1cs(F), [0]) and len(U.enumerate_witnesses(factor_r1cs(F), [0])) == 2 * 449 - 1


def test_roots_bac_hai():
    assert U._roots(449, 1, -1, 0, None) == {0, 1}          # u^2 - u = 0
    assert U._roots(449, 0, 0, 5, None) == set()
    assert U._roots(449, 0, 0, 0, None) is None
    assert U._roots(449, 0, 2, -6, None) == {3}
    assert U._roots(449, 1, 0, -4, {2, 5}) == {2}


def test_fuzz():
    c = mach_phan_tich()[0]
    r = c.compile_r1cs()
    x, w = c.r1cs_assignment({"N": 391, "a": 17, "b": 23})
    kq = U.fuzz_witness(r, x, w, trials=200, rng=random.Random(1))
    assert set(kq) == set(U.MUTATORS) and all(h == 0 for _, h, _ in kq.values())
    n, seen = U.fuzz_inputs(c, {"N": 391}, {"a": 17, "b": 23}, trials=200,
                            rng=random.Random(2024), spread={"a": (17, 23, 28, 30), "b": (17, 23, 28, 30)})
    assert set(seen) == {(23, 17), (28, 30), (30, 28)}


def test_z3_neu_co():
    if not U.HAVE_Z3:
        print("  (bỏ qua: không có z3)")
        return
    cu, moi = mach_phan_tich()[0], mach_phan_tich_sua()[0]
    for c, dung in ((cu, False), (moi, True)):
        r = c.compile_r1cs()
        x, w0 = c.r1cs_assignment({"N": 391, "a": 17, "b": 23})
        w2 = U.second_witness_z3(r, x, [int(v) for v in w0])
        assert (w2 is None) == dung
        if w2:
            assert tuple(w2[:2]) in {(23, 17), (28, 30), (30, 28)}
    assert U.all_witnesses_z3(cu.compile_r1cs(), [391]) == U.enumerate_witnesses(cu.compile_r1cs(), [391])


if __name__ == "__main__":
    _path.run_module(globals())
