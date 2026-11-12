"""Tỉ lệ lọt của giá trị mở sai trên trường toy F_449 (minh họa cho thm:ch3-batch-sound).

Mỗi lần thử: 4 đa thức bậc 15, tập điểm S = ([3,7], [3], [7,11], [3,7]); prover
trung thực tạo proof cho giá trị ĐÚNG, sau đó giá trị f_1(3) bị sửa thành +1 (kẻ gian
không thích nghi). Đếm số lần verifier vẫn chấp nhận. Với r = 449, kỳ vọng xấp xỉ
N/449: sai số tính đúng đắn của trường nhỏ là có thật, ở trường 255 bit thì ~2^-250.
Toy: Type 1 (G2 = G1).
"""
import os
import random
import sys
import time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from snarklib import kzg
from snarklib.poly import Poly

F = kzg.F
N = 2000
SRS = kzg.setup(15, 77)
Ss = [[3, 7], [3], [7, 11], [3, 7]]
rng = random.Random(7)
lot = {"batch-1d": 0, "gwc": 0, "shplonk": 0}
t0 = time.perf_counter()
for n in range(N):
    fs = [Poly(F, [rng.randrange(449) for _ in range(16)]) for _ in Ss]
    Cs = [kzg.commit(SRS, f) for f in fs]
    lab = b"thu-%d" % n
    ys, W = kzg.open_batch(SRS, fs, Cs, 3, kzg.fresh_transcript(lab))
    lot["batch-1d"] += kzg.verify_batch(SRS, Cs, 3, [ys[0] + 1] + ys[1:], W, kzg.fresh_transcript(lab))
    for ten, mo, kiem in (("gwc", kzg.open_gwc, kzg.verify_gwc),
                          ("shplonk", kzg.open_shplonk, kzg.verify_shplonk)):
        pr = mo(SRS, fs, Cs, Ss, kzg.fresh_transcript(lab))
        sai = [list(y) for y in pr.ys]
        sai[0][0] = sai[0][0] + 1
        lot[ten] += kiem(SRS, Cs, Ss, pr._replace(ys=sai), kzg.fresh_transcript(lab))
dt = time.perf_counter() - t0
print("N = %d lần thử, kỳ vọng N/449 = %.2f" % (N, N / 449))
for ten, c in lot.items():
    print("  %-9s lọt %d lần (%.3f%%)" % (ten, c, 100 * c / N))
print("thời gian: %.1f s" % dt)
assert all(c < 4 * N / 449 + 5 for c in lot.values())
