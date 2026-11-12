"""Chương 7 -- Sai số soundness của mini-PLONK + LogUp trên trường toy F_449.

Prover gian chứng minh (28, 30) với bảng {0..19}, {0..23} (bỏ phần dư của t).
Đo tỉ lệ proof gian được chấp nhận qua 200 nhãn transcript và tách phần do
ζ ∈ H (Z_H(ζ) = 0, mọi ràng buộc tại ζ 'tự thỏa'; xác suất 32/449 ≈ 7,1%).
Ở trường 256 bit phần này là n/|F|, không đáng kể.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib import logup

lc = logup.mach_phan_tich_lookup((19, 23))
pk, vk = logup.setup(lc, tau=123)
A, B, C, x = logup.assignment(391, 28, 30)
T = 200
lot = lot_H = zeta_H = 0
for i in range(T):
    lab = b"ch07/toy-%d" % i
    p = logup.LogupProver(pk, A, B, C, x, random.Random(i), label=lab, check=False)
    pf = p.run()
    inH = p.zeta ** p.n == 1
    zeta_H += inH
    if logup.verify(vk, x, pf, label=lab):
        lot, lot_H = lot + 1, lot_H + inH
print("proof gian được chấp nhận: %d/%d = %.1f%%" % (lot, T, 100 * lot / T))
print("  trong đó ζ ∈ H: %d; ζ ∉ H: %d" % (lot_H, lot - lot_H))
print("số lần ζ ∈ H (mọi proof): %d/%d = %.1f%%; lý thuyết 32/449 = %.1f%%"
      % (zeta_H, T, 100 * zeta_H / T, 100 * 32 / 449))
assert lot_H <= zeta_H and lot < T // 4
