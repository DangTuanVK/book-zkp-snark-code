"""Chương 12 -- Mạch đã sửa: thêm a <= 19, b <= 23 (lst:ch12-mach-da-sua).

Gadget assert_le (snarklib.dsl, thêm ở Ch.12): phân tích d = B - x thành 5 bit.
Kết hợp với a, b vừa 5 bit: a ∈ [0, 19], b ∈ [0, 23], nên a·b <= 437 < 449.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.dsl import mach_phan_tich, mach_phan_tich_sua
from snarklib import uniqueness as U

cu, moi = mach_phan_tich()[0], mach_phan_tich_sua(bounds=(19, 23))[0]
print("               R1CS m  |w|  Plonkish hàng  miền PLONK n   số witness (N = 391)")
kq = {}
for ten, c in (("mạch cũ (5 bit)", cu), ("mạch đã sửa", moi)):
    r1cs, pl = c.compile_r1cs(), c.compile_plonkish()
    ws = U.enumerate_witnesses(r1cs, [391])
    kq[ten] = (r1cs.m, pl.n, ws)
    print("%-17s %6d %4d %14d %13d   %d: %s" % (ten, r1cs.m, r1cs.n - 2, pl.n,
          1 << (pl.n - 1).bit_length(), len(ws), [tuple(w[:2]) for w in ws]))
print("Nhãn ràng buộc mạch đã sửa:", moi.stats()["labels"])

# Chứng minh bằng số: cận 19 = floor(sqrt(391)), 23 = floor(448/19)
import math
assert 19 == math.isqrt(391) and 23 == 448 // 19
print("max a·b = 19·23 = %d < 449 -> a·b ≡ 391 (mod 449) kéo theo a·b = 391 trên Z" % (19 * 23))
uoc = [(a, 391 // a) for a in range(2, 20) if 391 % a == 0 and 391 // a <= 23]
print("Ước của 391 với a <= 19, b <= 23:", uoc)

# Vì sao (23, 17) cũng biến mất: a <= 19 phá đối xứng giữa a và b
r_moi = moi.compile_r1cs()
for a, b in ((17, 23), (23, 17), (28, 30), (30, 28)):
    x, w = moi.r1cs_assignment({"N": 391, "a": a, "b": b})
    sai = [lab for (_, _, _, lab), r in zip(moi.constraints, r_moi.residuals(x, w)) if r]
    print("  (a, b) = (%2d, %2d): d_a = 19 - a ≡ %3d, d_b = 23 - b ≡ %3d (mod 449) -> %s"
          % (a, b, (19 - a) % 449, (23 - b) % 449,
             "thỏa" if not sai else "vi phạm %d ràng buộc %s" % (len(sai), sai)))
print("(23, 17) bị loại vì a = 23 > 19: cận không đối xứng chọn thứ tự a <= b;")
print("mệnh đề sửa đổi là 'a <= 19, b <= 23', khác mệnh đề gốc '2 <= a, b <= 31'.")

assert kq["mạch cũ (5 bit)"][0] == 13 and kq["mạch đã sửa"][0] == 25
assert [tuple(w[:2]) for w in kq["mạch đã sửa"][2]] == [(17, 23)] and uoc == [(17, 23)]
