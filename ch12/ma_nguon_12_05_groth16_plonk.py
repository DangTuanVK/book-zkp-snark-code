"""Chương 12 -- Groth16 và PLONK đều chấp nhận nghiệm giả (28, 30) trên mạch cũ.

Backend không có lỗi: proof đúng cho một witness thỏa mạch. Trên mạch đã sửa
Prover trung thực không sinh được witness, nên không có proof.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.curve import toy_generator
from snarklib.dsl import mach_phan_tich, mach_phan_tich_sua
from snarklib import groth16, plonk

E, P, G = toy_generator()                     # Type 1: G1 = G2 = G
gia = {"N": 391, "a": 28, "b": 30}            # 28·30 = 840 ≡ 391 (mod 449)
for ten, c in (("mạch cũ", mach_phan_tich()[0]), ("mạch đã sửa", mach_phan_tich_sua()[0])):
    r1cs = c.compile_r1cs()
    qap = groth16.r1cs_to_qap(r1cs)
    pk_g, vk_g, _ = groth16.setup(qap, G, G, random.Random(2016))
    circ = c.compile_plonkish()
    pk_p, vk_p = plonk.setup(circ, tau=123)
    print("== %s: R1CS m = %d (QAP n = %d), Plonkish %d hàng (PLONK n = %d) =="
          % (ten, r1cs.m, qap.n, circ.n, pk_p.n))
    for a, b in ((17, 23), (28, 30)):
        inp = {"N": 391, "a": a, "b": b}
        x, w = c.r1cs_assignment(inp)
        try:
            pf = groth16.prove(pk_g, qap, [1] + x + w, random.Random(a))
            kq_g = "proof -> verify = %s" % groth16.verify(vk_g, x, pf)
        except ValueError as e:
            kq_g = "Prover từ chối (%s)" % e
        try:
            pa, pb, pc, px = c.plonkish_assignment(inp)
            pf, _ = plonk.prove(pk_p, pa, pb, pc, px, random.Random(a))
            kq_p = "proof -> verify = %s" % plonk.verify(vk_p, px, pf)
        except ValueError as e:
            kq_p = "Prover từ chối (%s)" % e
        print("  (a, b) = (%d, %d)\n    Groth16: %s\n    PLONK  : %s" % (a, b, kq_g, kq_p))
        dung = c.compile_r1cs().is_satisfied(x, w)
        assert ("verify = True" in kq_g) == dung and ("verify = True" in kq_p) == dung
    assert c.compile_r1cs().is_satisfied(*c.r1cs_assignment(gia)) == (ten == "mạch cũ")
# Prover gian lận PLONK (check=False, bỏ phần dư) trên mạch đã sửa: 40 nhãn transcript
pa, pb, pc, px = c.plonkish_assignment(gia)
lot, zH = 0, 0
for i in range(40):
    lab = b"ch12/gian/%d" % i
    tr = {}
    pf, _ = plonk.prove(pk_p, pa, pb, pc, px, random.Random(i), label=lab, check=False)
    ok = plonk.verify(vk_p, px, pf, label=lab, trace=tr)
    lot += ok
    zH += tr["ZH"] == 0
print("Prover gian (28, 30), mạch đã sửa: lọt %d/40, trong đó ζ ∈ H (Z_H(ζ) = 0): %d;"
      " lý thuyết (hàng không vi phạm)/r = 62/449 ≈ %.1f%%" % (lot, zH, 6200 / 449))
print("Nghiệm giả là lỗi của MẠCH (thiếu ràng buộc), không phải của hệ chứng minh.")
