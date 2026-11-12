"""Chương 7 -- Phân tích bit và LogUp trên mạch N = 391 (tab:ch7-bit-vs-lookup).

Cột 1: mạch Ch.5-6 (a, b phân tích 5 bit) chạy mini-PLONK (snarklib.plonk).
Cột 2, 3: mạch 3 hàng + LogUp bảng {0..31} / hai bảng {0..19}, {0..23}.
Mọi chi phí đếm bằng snarklib.opcount; Prover chọn seed không phải dừng lại.
Phần cuối: M giá trị 5 bit cần kiểm phạm vi -> số hàng và cỡ miền n.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib import plonk, logup, opcount as oc
from snarklib.dsl import mach_phan_tich, Circuit

def do(tao, verify):
    """tao(seed) -> Prover; seed 2019 (như Ch.6) hoặc seed nhỏ nhất mà Prover không
    phải dừng lại và ζ not in H (ζ in H, xác suất 32/449, làm nhiều vô hướng bằng 0)."""
    def tot(sd):
        p = tao(sd); p.run()
        return p.aborts == 0 and p.zeta ** p.n != 1
    seed = next(sd for sd in [2019] + list(range(100)) if tot(sd))
    prover = tao(seed)
    with oc.scope() as cp:
        pf = prover.run()
    with oc.scope() as cv:
        assert verify(pf)
    assert prover.aborts == 0
    return cp, cv

c = mach_phan_tich()[0]
circ = c.compile_plonkish()
A, B, C, x = c.plonkish_assignment({"N": 391, "a": 17, "b": 23})
pk, vk = plonk.setup(circ, tau=123)
st = c.stats()
cot = [("phân tích 5 bit", st["m"], circ.n, pk.n, 8, "9 G1 + 6 F",
        do(lambda sd: plonk.Prover(pk, A, B, C, x, random.Random(sd)), lambda pf: plonk.verify(vk, x, pf)))]
for bounds, ten in [((31, 31), "LogUp {0..31}"), ((19, 23), "LogUp {0..19},{0..23}")]:
    lc = logup.mach_phan_tich_lookup(bounds)
    lpk, lvk = logup.setup(lc, tau=123)
    s = logup.circuit_stats(lc)
    A2, B2, C2, x2 = logup.assignment(391, 17, 23)
    pr = lambda sd: logup.LogupProver(lpk, A2, B2, C2, x2, random.Random(sd))
    cot.append((ten, "-", "%d + bảng %s" % (s["gate_rows"], "+".join(map(str, s["tables"]))),
                s["n"], s["preprocessed_comms"], "%d G1 + %d F" % (s["proof_G1"], s["proof_F"]),
                do(pr, lambda pf: logup.verify(lvk, x2, pf))))

dong = [("ràng buộc R1CS", lambda k: k[1]), ("hàng Plonkish", lambda k: k[2]),
        ("miền n", lambda k: k[3]), ("cam kết tiền xử lý", lambda k: k[4]),
        ("proof", lambda k: k[5]),
        ("Prover: số MSM / tổng cỡ", lambda k: "%d / %d" % (k[6][0]["msm"], k[6][0]["msm_size"])),
        ("Prover: nhân vô hướng", lambda k: k[6][0]["msm_nontrivial"]),
        ("Prover: NTT cỡ n", lambda k: k[6][0]["ntt"]),
        ("Verifier: MSM (cỡ)", lambda k: k[6][1]["msm_sizes"]),
        ("Verifier: nhân vô hướng", lambda k: k[6][1]["msm_nontrivial"]),
        ("Verifier: pairing", lambda k: k[6][1]["pairing"])]
print("%-26s" % "" + "".join("%-24s" % k[0] for k in cot))
for ten, fn in dong:
    print("%-26s" % ten + "".join("%-24s" % (fn(k),) for k in cot))

# M giá trị cần kiểm 5 bit: phân tích bit (DSL) và lookup (1 hàng mỗi giá trị + bảng 32)
print("\nM giá trị | bit: R1CS  hàng  n  | lookup: hàng  bảng  n   (* : n > 64, không có trên F_449)")
for M in (1, 2, 4, 8, 16, 32):
    d = Circuit()
    for i in range(M):
        d.to_bits(d.private("v%d" % i), 5, prefix="v%d_" % i)
    rows = d.compile_plonkish().n
    nb = plonk.domain_size(rows)
    print("%9d | %9d %5d %3d%s| %12d %5d %3d" % (M, d.stats()["m"], rows, nb, "*" if nb > 64 else " ",
                                                M, 32, plonk.domain_size(max(M, 32))))
assert cot[0][1] == 13 and cot[0][2] == 20 and cot[1][3] == 32
assert all(k[6][1]["pairing"] == 2 for k in cot)
