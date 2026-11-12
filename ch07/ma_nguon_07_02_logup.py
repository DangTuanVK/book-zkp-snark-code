"""Chương 7 -- LogUp đơn biến trên H (lst:ch7-logup, fig:ch7-logup).

Cột f = cột dây a của mạch có lookup (hàng 0: N, hàng 1: a, hàng 2: b),
selector q tra bảng ở hàng 1, 2; bảng t = {0..31} chiếm đúng n = 32 hàng; cột
bội m. Với thách thức λ, tổng tích lũy φ_{i+1} = φ_i + q_i/(λ+f_i) - m_i/(λ+t_i)
phải quay về 0 sau n bước (đồng nhất thức LogUp tại X = λ).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib import logup
from snarklib.logup import F
from snarklib.transcript import Transcript

n = 32
t = [F(j) for j in range(32)]                       # bảng {0..31}
q = [F(1) if i in (1, 2) else F(0) for i in range(n)]

def cot_va_boi(a, b):
    """Cột f (= cột a của mạch, đệm 0) và cột bội m (m_j = #{i: q_i = 1, f_i = t_j})."""
    f = [F(391), F(a), F(b)] + [F(0)] * (n - 3)
    m = [F(0)] * n
    for i in range(n):
        if q[i] == 1 and int(f[i]) < 32:
            m[int(f[i])] += 1                       # t_j = j nên vị trí của f_i là f_i
    return f, m

def tong_tich_luy(f, m, lam):
    phi = [F(0)]
    for i in range(n):
        buoc = -m[i] / (lam + t[i])
        if q[i] == 1:
            buoc += F(1) / (lam + f[i])
        phi.append(phi[-1] + buoc)
    return phi

def rang_buoc_C1(f, m, phi, lam, i):
    """(φ_{i+1} - φ_i)(λ + f_i)(λ + t_i) - q_i(λ + t_i) + m_i(λ + f_i) -- bằng 0 nếu bước đúng."""
    return ((phi[(i + 1) % n] - phi[i]) * (lam + f[i]) * (lam + t[i])
            - q[i] * (lam + t[i]) + m[i] * (lam + f[i]))

tr = Transcript(b"ch07/logup", 449)
for a, b in [(17, 23), (40, 21)]:                   # 40·21 = 840 ≡ 391, 40 ∉ bảng
    f, m = cot_va_boi(a, b)
    tr.absorb(b"f,m", f + m)
    lam, rut = tr.challenge(b"lambda"), 1
    while any(lam + v == 0 for v in t + f[1:3]):    # λ trùng cực điểm: rút lại
        lam, rut = tr.challenge(b"lambda"), rut + 1
    phi = tong_tich_luy(f, m, lam)
    print("(a, b) = (%d, %d), λ = %d (rút %d lần)" % (a, b, int(lam), rut))
    print("   i  q_i  f_i  t_i  m_i  φ_i")
    hang = sorted({0, 1, 2, 3, n - 1} | {j + d for j in range(n) if m[j] for d in (0, 1)})
    for i in hang:
        print("  %2d  %3d  %3d  %3d  %3d  %3d" % (i, int(q[i]), int(f[i]), int(t[i]),
                                                int(m[i]), int(phi[i])))
    print("  φ_n = φ_32 = %d -> quay về 0: %s" % (int(phi[n]), phi[n] == 0))
    vi_pham = [i for i in range(n) if rang_buoc_C1(f, m, phi[:n], lam, i) != 0]
    print("  hàng vi phạm C_1 khi φ là cột n giá trị (φ(ω^n) = φ(ω^0)):", vi_pham)
    ref = logup.running_sum(f, t, m, lam, q)          # đối chiếu snarklib.logup
    assert ref == phi
    assert (phi[n] == 0) == (a < 32 and b < 32) and (vi_pham == []) == (phi[n] == 0)
    assert logup.identity_at([f[1], f[2]], t, m, lam) == (phi[n] == 0)
