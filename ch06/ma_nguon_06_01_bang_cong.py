"""Chương 6 -- Bảng cổng PLONK của mạch N = 391 (tab:ch6-bang-cong) và hoán vị σ.

Mạch Plonkish 20 hàng do snarklib.dsl biên dịch, đệm tới n = 32 (ω_32 = 221),
coset H, k1·H, k2·H với k1 = 2, k2 = 3. Hàng 0 là cổng công khai (q_L = 1, PI_0 = -N).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib import plonk
from snarklib.dsl import mach_phan_tich

c = mach_phan_tich()[0]
circ = c.compile_plonkish()
a, b, cc, x = c.plonkish_assignment({"N": 391, "a": 17, "b": 23})
pk, vk = plonk.setup(circ, tau=123)
n, om, F = pk.n, pk.omega, plonk.F
print("hàng dùng = %d, n = %d, ω = %d, k1 = %d, k2 = %d" % (circ.n, n, int(om), int(pk.k1), int(pk.k2)))
sgn = lambda v: int(v) - 449 if int(v) > 224 else int(v)     # in hệ số âm cho dễ đọc
print(" j  qL  qR  qO  qM  qC  PI |   a    b    c")
PI = circ.public_input_column(x)
for j in range(circ.n):
    s = [sgn(pk.sel_vals[k][j]) for k in ("qL", "qR", "qO", "qM", "qC")]
    print("%2d %3d %3d %3d %3d %3d %4d | %4d %4d %4d" % (j, *s, sgn(PI[j]), a[j], b[j], cc[j]))
print("PI_0 = -391 ≡ 58 (mod 449); hàng 20..31: mọi selector = 0, dây = 0 (đệm)")
dem = {}
for j in range(circ.n):
    key = "cộng/tuyến tính" if pk.sel_vals["qM"][j] == 0 else "nhân"
    dem[key] = dem.get(key, 0) + 1
print("số cổng:", dem, "(hàng 0 là cổng công khai, tính vào nhóm tuyến tính)")

# hoán vị σ: các vòng (lớp ô bằng nhau) -- ô (cột, hàng) -> nhãn k·ω^j
ten = lambda u: "%s%d" % ("abc"[u // n], u % n)
print("ràng buộc sao chép: %d cặp -> %d vòng σ" % (len(circ.copies), len(pk.cycles)))
for cyc in pk.cycles:
    print("  (" + " ".join(ten(u) for u in cyc) + ")",
          " σ*:", ", ".join("%s->%d" % (ten(u), int(pk.sstar[u])) for u in cyc))
co_dinh = 3 * n - sum(len(cy) for cy in pk.cycles)
print("ô cố định của σ (không sao chép):", co_dinh)
assert circ.is_satisfied(a, b, cc, x) and n == 32 and int(om) == 221
assert dem["nhân"] == 11
# σ* là hoán vị của tập nhãn H' = H + k1H + k2H (3n nhãn phân biệt)
assert sorted(int(v) for v in pk.sstar) == sorted(int(v) for v in pk.ids)
assert len(set(int(v) for v in pk.ids)) == 3 * n
