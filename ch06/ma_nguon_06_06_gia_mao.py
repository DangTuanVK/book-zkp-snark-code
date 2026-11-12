"""Chương 6 -- Giả mạo kiểu Frozen Heart khi transcript KHÔNG băm đầu vào công khai.

Prover chỉ biết 17·23 = 391 nhưng xuất proof cho x* ≠ 391. Cơ chế (chi tiết ở
docstring của plonk.forge_frozen_heart): cộng Q(X) = q0 + q1·X vào t_lo, nhận ζ,
chọn x* = 391 - n(ζ - 1)·Q(ζ) để bù đúng phần lệch Z_H(ζ)·Q(ζ) qua PI(ζ).
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib import plonk
from snarklib.dsl import mach_phan_tich

c = mach_phan_tich()[0]
circ = c.compile_plonkish()
a, b, cc, x = c.plonkish_assignment({"N": 391, "a": 17, "b": 23})
pk, vk = plonk.setup(circ, tau=123)
# các N có witness trong mạch: a·b mod 449 với a, b ∈ [2, 31]
co_witness = {(u * v) % 449 for u in range(2, 32) for v in range(2, 32)}
print("số giá trị N (mod 449) có witness: %d / 449" % len(co_witness))

# (1) Một lần thử, transcript không băm x: proof hợp lệ cho một x* không chọn trước
pf, xs, lan = plonk.forge_frozen_heart(pk, a, b, cc, [391], None, random.Random(3))
ok_khong_bam = plonk.verify(vk, [xs], pf, hash_public=False)
ok_co_bam = plonk.verify(vk, [xs], pf, hash_public=True)
print("(1) một lần thử: x* = %d, có witness? %s" % (xs, xs in co_witness))
print("    verify không băm x: %s | verify có băm x: %s" % (ok_khong_bam, ok_co_bam))

# (2) Nhắm N = 389 (số nguyên tố, không có witness kể cả theo mod 449): mài (q0, q1)
pf2, xs2, lan2 = plonk.forge_frozen_heart(pk, a, b, cc, [391], 389, random.Random(3))
print("(2) nhắm 389, không băm x: %d lần thử, verify(389) = %s, verify(391) = %s"
      % (lan2, plonk.verify(vk, [389], pf2, hash_public=False),
         plonk.verify(vk, [391], pf2, hash_public=False)))

# (3) Đối chứng trên trường toy: mài cũng phá được bản CÓ băm x (sai số ~1/449/lần thử)
pf3, xs3, lan3 = plonk.forge_frozen_heart(pk, a, b, cc, [391], 389, random.Random(3),
                                          hash_public=True)
print("(3) nhắm 389, CÓ băm x: %d lần thử, verify(389) = %s  <- giới hạn của r = 449"
      % (lan3, plonk.verify(vk, [389], pf3)))

# (4) Thống kê số lần thử khi nhắm 389 (không băm), 100 seed làm mù
so_lan = [plonk.forge_frozen_heart(pk, a, b, cc, [391], 389, random.Random(100 + s))[2]
          for s in range(100)]
print("(4) 100 seed: số lần thử trung bình %.1f (kỳ vọng ≈ 449), min %d, max %d"
      % (sum(so_lan) / 100, min(so_lan), max(so_lan)))
assert lan == 1 and ok_khong_bam and not ok_co_bam and xs != 391
assert 389 not in co_witness and plonk.verify(vk, [389], pf2, hash_public=False)
