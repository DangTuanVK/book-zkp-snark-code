"""Chương 5 -- Tái ngẫu nhiên hóa proof Groth16 (tính dễ uốn yếu).

Từ π = (A, B, C) hợp lệ, ai cũng tạo được π' hợp lệ khác mà không cần witness:
    A' = r^{-1}·A,   B' = r·B + r·s·[δ]_2,   C' = C + s·A.
Kiểm: A'B' = AB + sδ·A và C'δ = Cδ + sA·δ, nên hai vế tăng cùng một lượng.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.curve import toy_generator
from snarklib.dsl import mach_phan_tich
from snarklib.groth16 import r1cs_to_qap, setup, prove, verify, rerandomize

E, P, G = toy_generator()
c = mach_phan_tich()[0]
qap = r1cs_to_qap(c.compile_r1cs())
x, w = c.r1cs_assignment({"N": 391, "a": 17, "b": 23})
rng = random.Random(2016)
pk, vk, tau = setup(qap, G, G, rng)
proof = prove(pk, qap, [1] + x + w, rng)       # cùng proof với 05_03

r, s = 5, 7                                    # do bất kỳ ai chọn
moi = rerandomize(vk, proof, r, s)
print("π  = A %s, B %s, C %s" % proof)
print("π' = A %s, B %s, C %s" % moi)
print("π' khác π:", moi != proof, "| verify(π'):", verify(vk, x, moi))

# Công thức "ngây thơ" chỉ đổi A, B (A' = r^{-1}A, B' = rB) thì vẫn đúng (s = 0);
# còn chỉ cộng vào C mà không bù B thì sai:
sai = proof._replace(C=proof.C + s * proof.A)
print("chỉ sửa C' = C + sA (không bù B):", verify(vk, x, sai))
assert verify(vk, x, moi) and moi != proof and not verify(vk, x, sai)
assert verify(vk, x, rerandomize(vk, proof, r, 0))
dem = sum(verify(vk, x, rerandomize(vk, proof, rr, ss))
          for rr, ss in [(2, 0), (3, 1), (448, 100), (17, 23)])
print("thêm 4 cặp (r, s) khác, số proof verify đúng:", dem)
assert dem == 4
