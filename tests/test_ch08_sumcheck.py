"""Test Chương 8: sum-check, zerocheck, productcheck, HyperPlonk tối giản."""
import random
import _path  # noqa: F401
from snarklib.field import GF
from snarklib.transcript import Transcript
from snarklib import sumcheck as sc, mle
from snarklib.dsl import mach_phan_tich

F = GF(449)
rng = random.Random(8)


def _rand_tables(K, l):
    return [[F(rng.randrange(449)) for _ in range(1 << l)] for _ in range(K)]


def test_interpolate_eval():
    p = lambda x: 3 * x ** 3 + 5 * x + 7
    evals = [F(p(i)) for i in range(4)]
    for x in (0, 2, 10, 448):
        assert sc.interpolate_eval(F, evals, x) == F(p(x))


def test_completeness_products():
    for l in range(1, 6):
        for K in range(1, 5):
            T = _rand_tables(K, l)
            g, d = sc.sum_of_products([(1, list(range(K)))])
            proof, r, fin = sc.prove(F, T, g, d, Transcript(b"c"))
            true = sum((g([t[i] for t in T]) for i in range(1 << l)), F(0))
            assert proof["claim"] == true
            assert fin == [mle.evaluate(F, t, r) for t in T]
            ok, rr = sc.verify_with_oracle(F, l, d, proof, Transcript(b"c"), g,
                                           lambda pt: [mle.evaluate(F, t, pt) for t in T])
            assert ok and rr == r


def test_general_expression_hypernova_style():
    # Σ_x eq(β,x)·(c0·(M0z)(M1z) + c1·(M2z)) + γ·(M3z)  — đa thức tùy ý của nhiều MLE
    l = 3
    T = _rand_tables(5, l)                       # bảng 0 = eq, 1..4 = M_j z
    g, d = sc.sum_of_products([(7, [0, 1, 2]), (-1, [0, 3]), (5, [4])])
    assert d == 3
    proof, r, fin = sc.prove(F, T, g, d, Transcript(b"hn"))
    ok, rr = sc.verify_with_oracle(F, l, d, proof, Transcript(b"hn"), g,
                                   lambda pt: [mle.evaluate(F, t, pt) for t in T])
    assert ok


def test_false_claim_and_tamper():
    T = _rand_tables(2, 3)
    g = lambda u: u[0] * u[1]
    proof, _, _ = sc.prove(F, T, g, 2, Transcript(b"x"))
    bad = dict(proof, claim=proof["claim"] + 1)
    assert not sc.verify(F, 3, 2, bad, Transcript(b"x"))[0]
    assert not sc.verify(F, 3, 2, proof, Transcript(b"x"), claim=proof["claim"] + 1)[0]
    # sửa một giá trị s_2(2): vòng kiểm hoặc kiểm cuối phải bác (nhiều nhãn, cho lọt ≤ 1)
    acc = 0
    for k in range(4):
        p, _, _ = sc.prove(F, T, g, 2, Transcript(b"t%d" % k))
        polys = [list(s) for s in p["polys"]]
        polys[1][2] = polys[1][2] + 1
        ok, rr = sc.verify_with_oracle(F, 3, 2, dict(p, polys=polys), Transcript(b"t%d" % k), g,
                                       lambda pt: [mle.evaluate(F, t, pt) for t in T])
        acc += ok
    assert acc <= 1
    # sai số bậc
    assert not sc.verify(F, 3, 3, proof, Transcript(b"x"))[0]


def test_mul_count_matches_formula():
    CF = sc.counting_field(F)
    for l, K, d in ((4, 4, 3), (5, 2, 2), (3, 3, 3)):
        T = [[CF(v) for v in t] for t in _rand_tables(K, l)]
        g, dd = sc.sum_of_products([(1, list(range(K)))])
        CF.counter["mul"] = 0
        sc.prove(CF, T, g, d, Transcript(b"m"), claim=CF(0))
        assert CF.counter["mul"] == sc.prover_cost(l, K, d, K - 1)


def test_zerocheck():
    l = 4
    a = _rand_tables(1, l)[0]
    b = _rand_tables(1, l)[0]
    c = [x * y for x, y in zip(a, b)]
    g = lambda u: u[0] * u[1] - u[2]
    proof, r, fin = sc.zerocheck_prove(F, [a, b, c], g, 2, Transcript(b"z"))
    ok, rr, final, eqv = sc.zerocheck_verify(F, l, 2, proof, Transcript(b"z"))
    assert ok and eqv * g([mle.evaluate(F, t, rr) for t in (a, b, c)]) == final
    c2 = list(c); c2[5] = c2[5] + 1
    acc = 0
    for k in range(4):
        p, _, _ = sc.zerocheck_prove(F, [a, b, c2], g, 2, Transcript(b"z%d" % k))
        ok, rr, final, eqv = sc.zerocheck_verify(F, l, 2, p, Transcript(b"z%d" % k))
        acc += bool(ok and eqv * g([mle.evaluate(F, t, rr) for t in (a, b, c2)]) == final)
    assert acc <= 1


def test_product_tree_and_productcheck():
    l = 3
    f = [F(rng.randrange(1, 449)) for _ in range(1 << l)]
    v = sc.product_tree(F, f)
    prod = F(1)
    for x in f:
        prod = prod * x
    assert v[(2 << l) - 2] == prod and v[(2 << l) - 1] == 0
    one = lambda u: F(1)
    ident = lambda u: u[0]
    proof, r, fin, vv = sc.productcheck_prove(F, [f], ident, one, 1, Transcript(b"p"))
    ok, s = sc.productcheck_verify(F, l, 1, proof, Transcript(b"p"), ident, one,
                                   lambda pt: mle.evaluate(F, vv, pt),
                                   lambda pt: [mle.evaluate(F, f, pt)])
    assert ok and s == prod
    # khai sai tích
    bad = dict(proof, s=prod + 1)
    ok2, _ = sc.productcheck_verify(F, l, 1, bad, Transcript(b"p"), ident, one,
                                    lambda pt: mle.evaluate(F, vv, pt),
                                    lambda pt: [mle.evaluate(F, f, pt)])
    assert not ok2


def test_hyperplonk_example():
    circ, _ = mach_phan_tich()
    pk = circ.compile_plonkish()
    a, b, c, x = circ.plonkish_assignment({"N": 391, "a": 17, "b": 23})
    T = sc.plonkish_tables(F, pk, a, b, c, x)
    proof, v = sc.hyperplonk_prove(F, T, Transcript(b"hp"))
    assert sc.hyperplonk_verify(F, T, proof, v, Transcript(b"hp"))
    # nghiệm giả (28, 30) thỏa mạch trên F_449 -> cũng qua
    a2, b2, c2, _ = circ.plonkish_assignment({"N": 391, "a": 28, "b": 30})
    T2 = sc.plonkish_tables(F, pk, a2, b2, c2, x)
    assert sc.hyperplonk_verify(F, T2, *sc.hyperplonk_prove(F, T2, Transcript(b"g")), Transcript(b"g"))
    # phá sao chép nhưng giữ cổng
    j = next(i for i in range(pk.n) if pk.sel["qM"][i] != 0 and pk.sel["qL"][i] == 0 and pk.sel["qC"][i] == 0)
    a3, c3 = list(a), list(c)
    a3[j] = F(a[j]) + 1
    c3[j] = -(pk.sel["qM"][j] * a3[j] * b[j]) / pk.sel["qO"][j]
    T3 = sc.plonkish_tables(F, pk, a3, b, c3, x)
    acc = sum(sc.hyperplonk_verify(F, T3, *sc.hyperplonk_prove(F, T3, Transcript(b"k%d" % k)),
                                   Transcript(b"k%d" % k)) for k in range(4))
    assert acc <= 1


if __name__ == "__main__":
    _path.run_module(globals())
