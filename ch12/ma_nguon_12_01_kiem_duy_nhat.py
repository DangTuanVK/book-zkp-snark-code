"""Chương 12 -- Kiểm witness duy nhất của mạch 5 bit Chương 5 (lst:ch12-kiem-duy-nhat).

Vét cạn có lan truyền trên TOÀN BỘ F_449^12 (không giả định a, b < 32): mỗi
ràng buộc còn một ẩn cho tập nghiệm hữu hạn, rẽ nhánh khi hết lan truyền.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.dsl import mach_phan_tich
from snarklib import uniqueness as U

c = mach_phan_tich()[0]                 # N công khai; a, b vừa 5 bit; a·b = N
r1cs = c.compile_r1cs()
ten = U.witness_names(c)                # w = (a, b, a0..a4, b0..b4)
print("R1CS: m = %d ràng buộc, |w| = %d biến; w = (%s)" % (r1cs.m, len(ten), ", ".join(ten)))

st = {}
ws = U.enumerate_witnesses(r1cs, [391], stats=st)
print("Vét cạn có lan truyền: %d nút tìm kiếm, rẽ nhánh trên %s"
      % (st["nodes"], [ten[j - 2] for j in st["branch_vars"]]))
print("Số witness thỏa với N = 391:", len(ws))
for w in ws:
    a, b = w[0], w[1]
    loai = "đáp án" if a * b == 391 else "NGHIỆM GIẢ (tràn trường)"
    print("  a = %2d = %s, b = %2d = %s | a·b = %3d = %d·449 + %3d  -> %s"
          % (a, "".join(map(str, w[6:1:-1])), b, "".join(map(str, w[11:6:-1])),
             a * b, a * b // 449, a * b % 449, loai))

# Câu hỏi "witness thứ hai?" (kiểu Picus): bỏ w0 = đáp án, còn witness nào khác?
w0 = [w for w in ws if w[:2] == [17, 23]][0]
khac = [w for w in ws if w != w0]
print("Có witness khác (17, 23)?", "CÓ: " + str([tuple(w[:2]) for w in khac]) if khac else "không")

# Kiểm chéo ở mức mạch: thử 32 x 32 cặp (a, b) với witness sinh bằng gợi ý
bf = U.brute_force_inputs(c, {"N": 391}, {"a": range(32), "b": range(32)}, r1cs)
print("Kiểm chéo vét cạn (a, b) ∈ [0, 31]^2: %d cặp thỏa" % len(bf))
assert [tuple(w[:2]) for w in ws] == [(17, 23), (23, 17), (28, 30), (30, 28)]
assert sorted((d["a"], d["b"]) for d in bf) == [tuple(w[:2]) for w in ws]
assert all(r1cs.is_satisfied([391], w) for w in ws) and st["complete"]
