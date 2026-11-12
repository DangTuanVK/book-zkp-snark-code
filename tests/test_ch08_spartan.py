"""Test Chương 8: mini-Spartan với PCS Ligero."""
import _path  # noqa: F401
from snarklib.field import GF
from snarklib.transcript import Transcript
from snarklib.dsl import mach_phan_tich, Circuit
from snarklib import spartan, ligero, sumcheck as sc

F = GF(449)


def _example():
    circ, _ = mach_phan_tich()
    r1 = circ.compile_r1cs()
    return circ, r1, spartan.pad_r1cs(r1)


def test_honest_and_fake_witness():
    circ, r1, inst = _example()
    assert (inst["size"], inst["s"]) == (16, 4)
    for ab, lab in (((17, 23), b"h"), ((28, 30), b"f"), ((23, 17), b"s")):
        x, w = circ.r1cs_assignment({"N": 391, "a": ab[0], "b": ab[1]})
        p = spartan.prove(inst, r1.z(x, w), Transcript(lab), pcs=ligero)
        assert spartan.verify(inst, [391], p, Transcript(lab), ligero)


def test_wrong_witness_and_public_input():
    circ, r1, inst = _example()
    x, w = circ.r1cs_assignment({"N": 391, "a": 17, "b": 24})
    acc = sum(spartan.verify(inst, [391], spartan.prove(inst, r1.z(x, w), Transcript(b"b%d" % k),
              pcs=ligero), Transcript(b"b%d" % k), ligero) for k in range(4))
    assert acc <= 1
    x, w = circ.r1cs_assignment({"N": 391, "a": 17, "b": 23})
    acc = 0
    for k in range(4):
        p = spartan.prove(inst, r1.z(x, w), Transcript(b"p%d" % k), pcs=ligero)
        acc += spartan.verify(inst, [392], p, Transcript(b"p%d" % k), ligero)   # đổi N
        q = dict(p, vW=F(p["vW"]) + 1)
        acc += spartan.verify(inst, [391], q, Transcript(b"p%d" % k), ligero)   # sửa W~(r_y)
        vA, vB, vC = p["vABC"]
        q = dict(p, vABC=(vA + 1, vB, vC))
        acc += spartan.verify(inst, [391], q, Transcript(b"p%d" % k), ligero)   # sửa v_A
    assert acc <= 1


def test_two_copies_size32():
    C = Circuit()
    vals = {}
    for i, (a_, b_) in enumerate(((17, 23), (23, 17))):
        N, a, b = C.public("N%d" % i), C.private("a%d" % i), C.private("b%d" % i)
        C.to_bits(a, 5, prefix="a%d_" % i)
        C.to_bits(b, 5, prefix="b%d_" % i)
        C.assert_mul(a, b, N)
        vals.update({"N%d" % i: 391, "a%d" % i: a_, "b%d" % i: b_})
    r1 = C.compile_r1cs()
    inst = spartan.pad_r1cs(r1)
    x, w = C.r1cs_assignment(vals)
    assert inst["size"] == 32 and r1.m == 26
    p = spartan.prove(inst, r1.z(x, w), Transcript(b"2"), pcs=ligero)
    assert spartan.verify(inst, [391, 391], p, Transcript(b"2"), ligero)


def test_proof_size_and_count():
    circ, r1, inst = _example()
    sz = spartan.proof_size(inst, ligero, t=8)
    assert (sz["rounds"], sz["queries"], sz["field_total"], sz["hash_total"]) == (8, 8, 72, 33)
    CF = sc.counting_field(F)
    x, w = circ.r1cs_assignment({"N": 391, "a": 17, "b": 23})
    tot, ph = spartan.prover_mul_count(CF, inst, [CF(v) for v in r1.z(x, w)], Transcript(b"c"))
    assert ph["sc1"] == sc.prover_cost(4, 4, 3, 2) and ph["sc2"] == sc.prover_cost(4, 2, 2, 1)
    assert tot == sum(ph.values())


if __name__ == "__main__":
    _path.run_module(globals())
