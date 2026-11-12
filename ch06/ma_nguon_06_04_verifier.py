"""Chương 6 -- Verifier PLONK 12 bước (lst:ch6-verifier) trên proof của N = 391.

Mã Verifier in trong sách: hàm verify của snarklib/plonk.py (bước 1-12 theo
gabizon2019plonk.txt dòng 1477-1529). File này in các giá trị trung gian và
chạy các phép thử âm. Trên F_449 một kiểm tra sai có thể lọt với xác suất
≈ 1/449, nên mỗi phép thử âm chạy với 4 nhãn transcript.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib import plonk
from snarklib.plonk import F
from snarklib.dsl import mach_phan_tich

c = mach_phan_tich()[0]
circ = c.compile_plonkish()
a, b, cc, x = c.plonkish_assignment({"N": 391, "a": 17, "b": 23})
pk, vk = plonk.setup(circ, tau=123)
proof, _ = plonk.prove(pk, a, b, cc, x, random.Random(2019))

tr = {}
ok = plonk.verify(vk, x, proof, trace=tr)
for k in ("beta", "gamma", "alpha", "zeta", "v", "u"):
    print("  bước 4  %-6s = %d" % (k, int(tr[k])))
print("  bước 5  Z_H(ζ) = ζ^32 - 1 = %d" % int(tr["ZH"]))
print("  bước 6  L_0(ζ) = %d" % int(tr["L0"]))
print("  bước 7  PI(ζ)  = -391·L_0(ζ) = %d" % int(tr["PI"]))
print("  bước 8  r_0    = %d" % int(tr["r0"]))
print("  bước 11 [E]_1  = %d·[1]_1" % int(tr["E"]))
print("  bước 12 hai pairing ->", ok)

nhan = [b"plonk/thu-1", b"plonk/thu-2", b"plonk/thu-3", b"plonk/thu-4"]
def dem_lot(tao):
    """Số nhãn (trên 4) mà proof sai vẫn được chấp nhận."""
    return sum(1 for lab in nhan if tao(lab))

# (1) sửa một dây: a ở hàng 19 (cổng nhân 17·23 = 391) thành 18; Prover gian bỏ phần dư
a_sai = list(a); a_sai[19] = F(18)
def sua_day(lab):
    pf, _ = plonk.prove(pk, a_sai, b, cc, x, random.Random(7), label=lab, check=False)
    return plonk.verify(vk, x, pf, label=lab)
# (2) đổi đầu vào công khai 391 -> 389 với proof trung thực
def doi_pi(lab):
    pf, _ = plonk.prove(pk, a, b, cc, x, random.Random(8), label=lab)
    return plonk.verify(vk, [389], pf, label=lab)
# (3) sửa một giá trị mở ā trong proof
def sua_mo(lab):
    pf, _ = plonk.prove(pk, a, b, cc, x, random.Random(9), label=lab)
    return plonk.verify(vk, x, pf._replace(a_bar=pf.a_bar + 1), label=lab)

kq = {"sửa dây a[19] = 18": dem_lot(sua_day), "đổi N = 389": dem_lot(doi_pi),
      "sửa ā + 1": dem_lot(sua_mo)}
for k, v in kq.items():
    print("  thử âm %-20s: lọt %d/4" % (k, v))
assert ok and all(v <= 1 for v in kq.values())
