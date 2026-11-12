"""Hiding KZG kiểu PolyCommit_Ped [KZG10 §3.3] trên đường cong toy (Mục 3.1, bài tập).

C = [f(τ)]_1 + [f̂(τ)]H với f̂ ngẫu nhiên. H sinh bằng băm-ra-điểm để không ai
chọn được log_G(H) — nhưng nhóm toy chỉ có 449 phần tử nên log này vét cạn
được ngay: tính che giấu ở đây chỉ để minh họa cấu trúc. Toy: Type 1.
"""
import os
import random
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from snarklib import kzg
from snarklib.poly import Poly

F = kzg.F
rng = random.Random(4)
srs = kzg.hiding_setup(3, 10)
H, so_lan_thu = kzg.hash_to_g1()
f = Poly(F, [23, 17, 0, 1])                       # cùng f của ví dụ 3.2a
C_tran = kzg.commit(srs, f)                        # KZG thường: tất định
fh1 = Poly(F, [rng.randrange(449) for _ in range(4)])
fh2 = Poly(F, [rng.randrange(449) for _ in range(4)])
C1, C2 = kzg.commit_hiding(srs, f, fh1), kzg.commit_hiding(srs, f, fh2)
print("H = %s (sau %d lần thử băm)" % (H, so_lan_thu))
print("KZG thường: C = %s (hai lần cam kết f cho cùng C)" % C_tran)
print("Hiding: C_1 = %s, C_2 = %s (cùng f, khác f̂)" % (C1, C2))

y, yh, W = kzg.open_hiding(srs, f, fh1, 1)
ok = kzg.verify_hiding(srs, C1, 1, y, yh, W)
print("mở tại z = 1: y = %d, ŷ = %d, W = %s, verify = %s" % (int(y), int(yh), W, ok))
print("verify với y sai lệch 1 =", kzg.verify_hiding(srs, C1, 1, y + 1, yh, W))

beta = kzg.dlog_bruteforce(srs.G1, H)
print("HẠN CHẾ TOY: vét cạn 449 phần tử tìm được log_G(H) = %d" % beta)
# biết β, kẻ gian mở C_1 thành giá trị y' tùy ý: chọn ŷ' = ŷ + (y - y')/β
y_gia = F(100)
yh_gia = yh + (y - y_gia) / F(beta)
print("biết β: mở C_1 tại z = 1 thành y' = 100 với ŷ' = %d -> verify = %s"
      % (int(yh_gia), kzg.verify_hiding(srs, C1, 1, y_gia, yh_gia, W)))
assert C1 != C2 and ok and beta * srs.G1 == H
assert kzg.verify_hiding(srs, C1, 1, y_gia, yh_gia, W)
