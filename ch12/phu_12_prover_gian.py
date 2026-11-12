"""Chương 12 (phụ, không in trong sách) -- Vì sao Prover gian PLONK lọt 7/40 lần.

Lặp lại đúng 40 lần thử của ma_nguon_12_05 (mạch đã sửa, (a, b) = (28, 30),
check=False). Prover gian bỏ phần dư rem(X) của phép chia tử số cho Z_H; phép
kiểm của Verifier tại ζ đúng khi và chỉ khi rem(ζ) = 0 hoặc Z_H(ζ) = 0. Ta đếm
nghiệm của rem trong H và ngoài H, và đối chiếu với các lần lọt.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.dsl import mach_phan_tich_sua
from snarklib import plonk

c = mach_phan_tich_sua()[0]
pk, vk = plonk.setup(c.compile_plonkish(), tau=123)
pa, pb, pc, px = c.plonkish_assignment({"N": 391, "a": 28, "b": 30})
H = {pow(int(pk.omega), j, 449) for j in range(pk.n)}
lot, trong_H, nghiem_ngoai, cau_truc = [], [], [], set()
for i in range(40):
    lab = b"ch12/gian/%d" % i
    p = plonk.Prover(pk, pa, pb, pc, px, random.Random(i), lab, True, False)
    pf = p.run()
    _, rem = p.quotient_naive()
    R = {t for t in range(449) if int(rem(t)) == 0}
    cau_truc.add((rem.degree(), len(R & H), len(R - H)))
    z = int(p.zeta)
    ok = plonk.verify(vk, px, pf, label=lab)
    assert ok == (z in R or z in H)
    if ok:
        lot.append(i)
        (trong_H if z in H else nghiem_ngoai).append((i, z))
print("n = %d; (bậc rem, số nghiệm trong H, số nghiệm ngoài H) qua 40 lần: %s" % (pk.n, sorted(cau_truc)))
print("lọt %d/40: ζ ∈ H ở %d lần %s; ζ ∉ H nhưng rem(ζ) = 0 ở %d lần %s"
      % (len(lot), len(trong_H), [i for i, _ in trong_H], len(nghiem_ngoai), nghiem_ngoai))
print("xác suất lọt mỗi lần = (|H| + 1)/r = 65/449 ≈ %.1f%%" % (6500 / 449))
