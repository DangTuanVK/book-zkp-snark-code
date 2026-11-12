"""Chương 5 -- Số ràng buộc của mạch N = 391 ở ba dạng R1CS, Plonkish, CCS.

Dữ liệu cho bảng tab:ch5-so-rang-buoc. Bản "tường minh" thêm hai cổng
a - 1 != 0, b - 1 != 0 (gợi ý nghịch đảo) để so sánh chi phí.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.dsl import mach_phan_tich
from snarklib.groth16 import r1cs_to_qap

dung = {"N": 391, "a": 17, "b": 23}
for ten, nz in (("gọn (mặc định)", False), ("tường minh a,b >= 2", True)):
    c = mach_phan_tich(nonzero=nz)[0]
    st = c.stats()
    r1cs = c.compile_r1cs()
    x, w = c.r1cs_assignment(dung)
    pk = c.compile_plonkish()
    pa, pb, pc, px = c.plonkish_assignment(dung)
    ccs = c.compile_ccs()
    qap = r1cs_to_qap(r1cs)
    so_nhan = sum(1 for v in pk.sel["qM"] if v != 0)
    print("== Mạch", ten, "==")
    print("R1CS    : m = %d, n = %d (ℓ = %d, witness = %d), bậc 2, nnz(A,B,C) = %s, %s"
          % (st["m"], st["n"], st["n_public"], st["n_witness"], st["nnz"], st["labels"]))
    print("Plonkish: %d hàng (đệm 2^k: %d), %d cổng nhân, %d cổng tuyến tính/công khai,"
          " %d dây tạm, %d ràng buộc sao chép, bậc cổng 2"
          % (pk.n, 1 << (pk.n - 1).bit_length(), so_nhan, pk.n - so_nhan, pk.n_temps, len(pk.copies)))
    print("CCS     : m = %d, n = %d, t = %d ma trận, q = %d đa tập %s, c = %s, bậc d = %d"
          % (ccs.m, ccs.n, len(ccs.M), len(ccs.S), ccs.S, [int(v) if int(v) < 225 else int(v) - 449 for v in ccs.c], ccs.degree()))
    print("QAP/Groth16: m + ℓ + 1 = %d hàng -> miền H cỡ n = %d" % (st["m"] + st["n_public"] + 1, qap.n))
    ok = (r1cs.is_satisfied(x, w), pk.is_satisfied(pa, pb, pc, px), ccs.is_satisfied(x, w))
    print("Witness (17, 23) thỏa R1CS / Plonkish / CCS:", ok)
    assert all(ok)

c0 = mach_phan_tich()[0]
assert c0.stats()["m"] == 13 and c0.compile_plonkish().n == 20
# a, b >= 2 là hệ quả: a = 0 -> a·b = 0 != 391; a = 1 -> b ≡ 391 không vừa 5 bit
assert all((1 * b) % 449 != 391 for b in range(32))
print("Kiểm: với a, b in [0, 31], a in {0, 1} không thể cho a·b ≡ 391 (mod 449)")
