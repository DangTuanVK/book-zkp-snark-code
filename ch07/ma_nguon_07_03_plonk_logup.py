"""Chương 7 -- mini-PLONK + LogUp: mạch a·b = N kiểm phạm vi bằng bảng.

Mạch 3 hàng (công khai N, cổng nhân, hàng chép b sang cột a) + lookup cột a.
bounds = (31, 31): một bảng {0..31} (thay phân tích 5 bit của Ch.5-6).
bounds = (19, 23): hai bảng {0..19}, {0..23} (bản sửa của Ch.12).
Thử âm chạy 12 nhãn transcript: trên F_449 proof gian lọt khi ζ rơi vào H
(xác suất 32/449 ≈ 7%) hoặc trùng hợp ≈ 1/449 khác.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib import logup
from snarklib.logup import F

NHAN = [b"ch07/thu-%d" % i for i in range(12)]

def chay(bounds, a, b, check=True, seed=2, label=logup.LABEL):
    lc = logup.mach_phan_tich_lookup(bounds)
    pk, vk = logup.setup(lc, tau=123)
    A, B, C, x = logup.assignment(391, a, b)
    p = logup.LogupProver(pk, A, B, C, x, random.Random(seed), label=label, check=check)
    pf = p.run()
    return logup.verify(vk, x, pf, label=label), pf, p, vk

for bounds in [(31, 31), (19, 23)]:
    ok, pf, p, vk = chay(bounds, 17, 23)
    g1, fr = logup.proof_size(pf)
    print("bounds = %s: n = %d, K = %d lookup, bảng %s" % (bounds, p.n, len(vk.lk),
          [len(L["table_raw"]) for L in p.pk.lk]))
    print("  (17, 23): verify = %s; proof %d G1 + %d F; λ = %d, ζ = %d; dừng lại %d lần"
          % (ok, g1, fr, int(p.lam), int(p.zeta), p.aborts))
    print("  bội m khác 0:", [{int(L["table"][j]): int(v) for j, v in enumerate(m) if v}
                              for L, m in zip(p.pk.lk, p.m_vals)])
    assert ok and p.zeta ** p.n != 1

# nghiệm giả (28, 30): bảng {0..31} vẫn nhận (giống phân tích 5 bit), bảng (19, 23) bác
ok31 = chay((31, 31), 28, 30)[0]
print("(28, 30) với bảng {0..31}: verify =", ok31, "(28·30 = 840 ≡ 391)")
try:
    chay((19, 23), 28, 30)
except ValueError as e:
    print("(28, 30) với bảng (19, 23): Prover trung thực dừng --", e)

def dem_lot(bounds, a, b, sua=None):
    """(số nhãn lọt, trong đó bao nhiêu lần ζ in H, tức ζ^n = 1)."""
    lot = lot_H = 0
    for i, lab in enumerate(NHAN):
        _, pf, p, vk = chay(bounds, a, b, check=sua is not None, seed=50 + i, label=lab)
        if sua is not None:
            pf = sua(pf)
        if logup.verify(vk, [391], pf, label=lab):
            lot, lot_H = lot + 1, lot_H + (p.zeta ** p.n == 1)
    return lot, lot_H

sua_m = lambda pf: pf._replace(lk=[pf.lk[0]._replace(m_bar=pf.lk[0].m_bar + 1)] + pf.lk[1:])
kq = {"(28, 30), bảng (19, 23), Prover gian": dem_lot((19, 23), 28, 30),
      "(40, 21), bảng {0..31}, Prover gian": dem_lot((31, 31), 40, 21),
      "(17, 23), sửa m̄ + 1": dem_lot((31, 31), 17, 23, sua_m)}
for k, v in kq.items():
    print("  thử âm %-38s: lọt %d/%d (ζ in H: %d)" % (k, v[0], len(NHAN), v[1]))
assert ok31 and all(v[0] <= 3 for v in kq.values())
