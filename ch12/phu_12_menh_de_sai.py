"""Chương 12 (phụ, không in trong sách) -- Mệnh đề SAI chứng minh được trên mạch cũ.

Với mạch 5 bit của Chương 5, đếm các N ∈ [1, 448] mà mọi witness đều không
phải phân tích thật (a·b = N trên Z với 2 <= a, b): khi đó mệnh đề "N có hai
thừa số trong [2, 31]" sai nhưng vẫn có proof hợp lệ. Vét cạn đầy đủ bằng
uniqueness.enumerate_witnesses trên toàn F_449^12.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.dsl import mach_phan_tich
from snarklib import uniqueness as U

r1cs = mach_phan_tich()[0].compile_r1cs()
sai, gia = [], []
for N in range(1, 449):
    ab = [tuple(w[:2]) for w in U.enumerate_witnesses(r1cs, [N])]
    that = [(a, b) for a, b in ab if a * b == N and a >= 2 and b >= 2]
    if any(a * b >= 449 for a, b in ab):
        gia.append(N)
        if not that:
            sai.append((N, ab))
print("Số N có nghiệm giả:", len(gia), "/ 448")
print("Số N có nghiệm giả nhưng KHÔNG có phân tích thật trong [2, 31]^2:", len(sai))
for N, ab in sai:
    if N in (37, 47):
        print("  N = %d (số nguyên tố): witness %s" % (N, ab))
assert len(gia) == 93
