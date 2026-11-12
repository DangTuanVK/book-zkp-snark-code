"""Chương 12 -- Bỏ/thêm từng ràng buộc phạm vi: số witness của N = 391 thay đổi ra sao.

Mỗi biến thể được vét cạn có lan truyền trên toàn F_449^|w| (đầy đủ), so với
bản lookup của Chương 7 (snarklib.logup, cận theo bảng).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.dsl import Circuit
from snarklib import uniqueness as U, logup


def mach(bit_a, bit_b, le_a, le_b):
    c = Circuit()
    N, a, b = c.public("N"), c.private("a"), c.private("b")
    if bit_a: c.to_bits(a, 5, prefix="a")
    if bit_b: c.to_bits(b, 5, prefix="b")
    if le_a: c.assert_le(a, 19, prefix="da")
    if le_b: c.assert_le(b, 23, prefix="db")
    c.assert_mul(a, b, N)
    return c


BIEN_THE = [("chỉ a·b = N", (0, 0, 0, 0)),
            ("5 bit a, b (mạch cũ Ch.5)", (1, 1, 0, 0)),
            ("chỉ a <= 19, b <= 23 (không 5 bit)", (0, 0, 1, 1)),
            ("5 bit a, b + a <= 19", (1, 1, 1, 0)),
            ("5 bit a, b + a <= 19, b <= 23 (đã sửa)", (1, 1, 1, 1))]
print("%-40s %4s %6s  %s" % ("biến thể", "m", "hàng", "số witness / vài (a, b)"))
kq = {}
for ten, cfg in BIEN_THE:
    c = mach(*cfg)
    r1cs = c.compile_r1cs()
    ws = U.enumerate_witnesses(r1cs, [391])
    ab = [tuple(w[:2]) for w in ws]
    kq[cfg] = ab
    print("%-40s %4d %6d  %3d  %s%s" % (ten, r1cs.m, c.compile_plonkish().n, len(ab),
          ab[:6], " ..." if len(ab) > 6 else ""))
# Duy nhất với N = 391 chưa đủ: đếm mọi N ∈ [1, 448] có nghiệm giả (a·b >= 449 trên Z)
def mien(bit, le, B):
    d = set(range(32)) if bit else set(range(449))
    return d & {(B - k) % 449 for k in range(32)} if le else d


def so_N_gia(cfg):
    A, Bm = mien(cfg[0], cfg[2], 19), mien(cfg[1], cfg[3], 23)
    return len({a * b % 449 for a in A for b in Bm if a * b >= 449 and a * b % 449})


for ten, cfg in BIEN_THE:
    print("  %-40s số N có nghiệm giả: %3d / 448" % (ten, so_N_gia(cfg)))
c20 = mach(0, 0, 1, 1)
ws20 = [tuple(w[:2]) for w in U.enumerate_witnesses(c20.compile_r1cs(), [20])]
print("Ví dụ không 5 bit, N = 20: %s; 445 = -4, 444 = -5 (mod 449)" % ws20)
print("Bản lookup Ch.7: bảng {0..31}: %s; bảng {0..19}, {0..23}: %s"
      % (logup.vet_can(391, (31, 31)), logup.vet_can(391, (19, 23))))
st = logup.circuit_stats(logup.mach_phan_tich_lookup((19, 23)))
print("  lookup (19, 23): %d hàng cổng, %d tra bảng, bảng %s, n = %d"
      % (st["gate_rows"], st["lookups"], st["tables"], st["n"]))
assert len(kq[(0, 0, 0, 0)]) == 448 and len(kq[(1, 1, 0, 0)]) == 4
assert kq[(1, 1, 1, 1)] == [(17, 23)] and so_N_gia((1, 1, 1, 1)) == 0
assert (445, 444) in ws20 and logup.vet_can(391, (19, 23)) == [(17, 23)]
