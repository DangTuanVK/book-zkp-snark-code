"""Chương 12 -- Hỏi z3: "tồn tại witness thứ hai?" (kiểu Picus) cho hai mạch.

Chạy trong môi trường có z3-solver (env zksnark, xem ket_qua/ch12/z3_env.txt).
Ràng buộc trên F_449 được mã hóa bằng bit-vector 32 bit với phép URem 449.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.dsl import mach_phan_tich, mach_phan_tich_sua
from snarklib import uniqueness as U

if not U.HAVE_Z3:
    sys.exit("không có z3: chạy ma_nguon_12_01 (vét cạn thuần Python) thay thế")
print("z3", U.z3.get_version_string())
for ten, c in (("mạch cũ", mach_phan_tich()[0]), ("mạch đã sửa", mach_phan_tich_sua()[0])):
    r1cs = c.compile_r1cs()
    x, w0 = c.r1cs_assignment({"N": 391, "a": 17, "b": 23})
    w0 = [int(v) for v in w0]
    t = time.time()
    w2 = U.second_witness_z3(r1cs, x, w0)          # có w != w0 cùng x không?
    dt = time.time() - t
    if w2 is None:
        print("%-11s: unsat sau %.1f s -> witness của N = 391 là DUY NHẤT" % (ten, dt))
    else:
        print("%-11s: sat sau %.1f s -> witness thứ hai (a, b) = (%d, %d)" % (ten, dt, w2[0], w2[1]))
    t = time.time()
    tat_ca = U.all_witnesses_z3(r1cs, x)           # chặn dần tới khi unsat
    print("             liệt kê bằng mệnh đề chặn: %d witness %s (%.1f s)"
          % (len(tat_ca), [tuple(w[:2]) for w in tat_ca], time.time() - t))
    assert tat_ca == U.enumerate_witnesses(r1cs, x)  # khớp vét cạn thuần Python
    assert (w2 is None) == (ten == "mạch đã sửa")
