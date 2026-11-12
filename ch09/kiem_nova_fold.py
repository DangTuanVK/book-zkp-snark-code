"""Chương 9 — Kiểm folding Nova cho relaxed R1CS trên mạch ví dụ (CHỈ lấy dữ liệu).

Công thức gộp đã có ở [ZKP, Mục 13.1] (kothapalli2022nova.txt:745-748, :798-806);
sách không in lại mã này, chỉ dùng số liệu: một bước gộp, gộp hai instance đã
relaxed (E ≠ 0), kiểm âm. Không cam kết: T, E, W gửi rõ.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.dsl import mach_phan_tich
from snarklib.folding import RelaxedR1CS
from snarklib.transcript import Transcript

c = mach_phan_tich()[0]
R = RelaxedR1CS(c.compile_r1cs())
wit = lambda a, b: R.trivial(*c.r1cs_assignment({"N": 391, "a": a, "b": b}))
(i1, W1), (i2, W2), (i3, W3) = wit(17, 23), wit(23, 17), wit(17, 23)

I, W, T, r = R.fold(i1, W1, i2, W2, Transcript(b"ch9-nova"))
print("Bước 1: r =", r, " u =", I["u"], " x =", I["x"], " thỏa:", R.is_satisfied(I, W))
print("  T =", T)
print("  E =", I["E"], "(%d phần tử khác 0)" % sum(1 for v in I["E"] if v))
assert R.is_satisfied(I, W) and I["u"] == 1 + r
I2, W2n, T2, r2 = R.fold(I, W, i3, W3, Transcript(b"ch9-nova-2"))
print("Bước 2 (relaxed + mới): r =", r2, " u =", I2["u"], " thỏa:", R.is_satisfied(I2, W2n))
assert R.is_satisfied(I2, W2n)
print("Kích thước: T có m = %d phần tử (bản thật: 1 cam kết); instance (E, u, x): %d + 1 + 1"
      % (len(T), len(I["E"])))

def thu(a, b, k=20):
    ib, Wb = wit(a, b)
    return sum(R.is_satisfied(*R.fold(i1, W1, ib, Wb, Transcript(b"neg%d" % i))[:2]) for i in range(k))
sai, gia = thu(17, 24), thu(28, 30)
print("Gộp với (17, 24): thỏa %d/20 | gộp với nghiệm giả (28, 30): thỏa %d/20" % (sai, gia))
assert sai <= 2 and gia == 20
