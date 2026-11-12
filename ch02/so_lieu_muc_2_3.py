"""Chương 2, Mục 2.3 -- số lần gọi hàm băm trung bình của ba tấn công trong
Mã nguồn 2.2 (tệp phụ, không in): FS yếu không nhắm, FS yếu nhắm N* = 43,
FS mạnh với N = 43. So với kỳ vọng 449/224, 449 và 449/4."""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ma_nguon_02_02_weak_fs import (F, Poly, X1, C6, gia_mao_yeu, gia_mao_manh,
                                    co_witness, r)

rng = random.Random(7)
khong_wit = sum(1 for N in range(r) if not co_witness(N))
print("số N không có witness trong R6: %d / %d" % (khong_wit, r))
LAP = 200
yeu, nham, manh = [], [], []
for _ in range(LAP):
    A = Poly(F, [17]) + Poly(F, [rng.randrange(r)]) * X1
    Q = (A * (A + C6) - Poly(F, [391])) // X1
    yeu.append(gia_mao_yeu(A, Q, None)[2])
    nham.append(gia_mao_yeu(A, Q, 43)[2])
    manh.append(gia_mao_manh(43, rng)[2])
tb = lambda v: sum(v) / len(v)
print("%d lần lặp lại mỗi tấn công:" % LAP)
print("  FS yếu, không nhắm: trung bình %.2f lần gọi băm (kỳ vọng %.2f)"
      % (tb(yeu), r / khong_wit))
print("  FS yếu, nhắm N* = 43: trung bình %.1f (kỳ vọng %d)" % (tb(nham), r))
print("  FS mạnh, N = 43: trung bình %.1f (kỳ vọng %.2f), nhỏ nhất %d, lớn nhất %d"
      % (tb(manh), r / 4, min(manh), max(manh)))
