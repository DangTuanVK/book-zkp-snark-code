"""Chương 7 -- Kiểm đa tập bằng tích (eq:ch7-multiset) trên F_449.

Tra bảng f = (a, b) = (17, 23) vào t = {0..31}: bội m_17 = m_23 = 1. Kiểm
Π (X + f_i) = Π (X + t_j)^{m_j} như đa thức và tại thách thức γ ngẫu nhiên;
đổi f thành (17, 40) (40 ∉ t) thì đẳng thức sai. Tích tích lũy Z_i là cơ chế
z(X) của PLONK (Ch.6) áp cho hai dãy f và t^m.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib import multiset as ms
from snarklib.multiset import F449 as F

t = list(range(32))                     # bảng kiểm phạm vi 5 bit
f = [17, 23]                            # a, b của ví dụ xuyên suốt
m, ngoai = ms.multiplicities(f, t)
print("f =", f, " bội khác 0:", {t[j]: k for j, k in enumerate(m) if k}, " ngoài bảng:", ngoai)

# (1) đồng nhất thức đa thức: vế phải bỏ các t_j có m_j = 0
trai = ms.product_poly(F, f)
phai = ms.product_poly(F, t, m)
print("Π(X + f_i)        =", trai)
print("Π(X + t_j)^{m_j}  =", phai, " -> bằng nhau:", trai == phai)

# (2) Verifier chỉ kiểm tại γ ngẫu nhiên (Schwartz–Zippel: lọt ≤ deg/p)
rng = random.Random(7)
gamma = rng.randrange(449)
print("γ = %d: hai vế = %d, %d" % (gamma, int(ms.product_at(F, f, gamma)),
                                  int(ms.product_at(F, t, gamma, m))))

# (3) giá trị ngoài bảng: f' = (17, 40), Prover gian vẫn khai bội m như trên
f2 = [17, 40]
g = ms.expand(t, m)                     # t^m = (17, 23)
_, ngoai2 = ms.multiplicities(f2, t)
hits, tong = ms.false_positive_rate(F, f2, g)
print("f' = %s: ngoài bảng %s; số γ ∈ F_449 làm hai tích bằng nhau: %d/%d (cận deg/p = %d/449)"
      % (f2, ngoai2, hits, tong, len(f2)))

# (4) tích tích lũy Z_0 = 1, Z_{i+1} = Z_i (γ + f_i)/(γ + g_i): không phụ thuộc thứ tự
Z = ms.grand_product(F, f, g[::-1], gamma)
Z2 = ms.grand_product(F, f2, g, gamma)
print("f = (17, 23), g = (23, 17): Z =", [int(z) for z in Z], "-> Z_n = 1:", Z[-1] == 1)
print("f' = (17, 40), g = (17, 23): Z =", [int(z) for z in Z2], "-> Z_n = 1:", Z2[-1] == 1)

# (5) vì sao m là số mũ thì khó: bậc vế phải = Σ m_j = n, không phải bậc theo hàng
print("bậc Π(X + t_j)^{m_j} = Σ m_j =", phai.degree())

assert trai == phai and ngoai == [] and ngoai2 == [40]
assert ms.lookup_product_check(F, f, t, m, gamma) and Z[-1] == 1
assert not ms.multiset_equal_poly(F, f2, g) and hits <= len(f2) and Z2[-1] != 1
