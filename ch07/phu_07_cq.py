"""Chương 7 -- Kiểm chứng công thức thương lưu sẵn (eq:ch7-cq) trên F_449.

Bảng t = {0..31} đặt trên miền V = <ω_32>, hai phép tra f = (17, 23) đặt trên
miền H = <ω_2> = {1, -1}. Với thách thức β (rút từ transcript, tránh cực điểm):
  A_i = m_i/(β + t_i) trên V,   B_j = 1/(β + f_j) trên H.
Kiểm: (1) Σ_V A = Σ_H B, đọc qua hệ số tự do: |V|·A(0) = |H|·B(0);
(2) A(X)(T(X) + β) - m(X) chia hết cho Z_V(X), thương Q_A;
(3) Q_A = Σ_{i: A_i ≠ 0} A_i·Q_i với Q_i = L_i(X)(T(X) - t_i)/Z_V(X) tính sẵn,
    nên [Q_A]_1 là MSM cỡ = số bội khác 0 thay vì cỡ |V| - 1.
(4) Caulk: với tập chỉ số I = {17, 23}, T(X) - C_I(X) = Z_I(X)·H_I(X) và
    H_I = Σ_{i in I} μ_i·W_i, μ_i = 1/Z_I'(ω^i), W_i = (T(X) - t_i)/(X - ω^i)
    là bằng chứng mở KZG tính sẵn tại ω^i.
(5) Lasso: bảng {0..1023} phân rã thành hai bảng con {0..31}:
    T[(r1, r2)] = 32·T1[r1] + T2[r2]; MLE của bảng con đồng nhất (5 biến,
    big-endian) là Σ 2^(5-s) x_s.
Đây là kiểm chứng đại số của các công thức trong sách, không phải cài đặt
giao thức cq, Caulk hay Lasso.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib import kzg, opcount as oc
from snarklib.field import GF
from snarklib.poly import Poly, lagrange_interpolate, vanishing_poly, domain, lagrange_basis
from snarklib.transcript import Transcript

F = GF(449)
NV, NH = 32, 2
V, H = domain(F, NV), domain(F, NH)
t = [F(i) for i in range(NV)]
f = [F(17), F(23)]
m = [F(sum(1 for v in f if v == ti)) for ti in t]

tr = Transcript(b"ch07/cq", 449)
tr.absorb(b"f,m", f + m)
beta, rut = tr.challenge(b"beta"), 1
while any(beta + v == 0 for v in t + f):
    beta, rut = tr.challenge(b"beta"), rut + 1
print("|V| = %d, |H| = %d, beta = %d (rút %d lần)" % (NV, NH, int(beta), rut))

A_vals = [mi / (beta + ti) for mi, ti in zip(m, t)]
B_vals = [F(1) / (beta + fj) for fj in f]
A = lagrange_interpolate(F, V, A_vals)
B = lagrange_interpolate(F, H, B_vals)
T = lagrange_interpolate(F, V, t)
M = lagrange_interpolate(F, V, m)
nz = [i for i in range(NV) if A_vals[i] != 0]
print("A_i khác 0 tại i =", nz, " giá trị:", [int(A_vals[i]) for i in nz])
print("B_j =", [int(b) for b in B_vals])
print("tổng trên V = %d, tổng trên H = %d" % (int(sum(A_vals, F(0))), int(sum(B_vals, F(0)))))
print("|V|·A(0) = %d, |H|·B(0) = %d" % (int(F(NV) * A(F(0))), int(F(NH) * B(F(0)))))

ZV = vanishing_poly(F, NV)
QA, du = (A * (T + beta) - M).divmod(ZV)
print("A(X)(T(X)+beta) - m(X) chia hết cho Z_V:", du.is_zero(), "; bậc Q_A =", QA.degree())

Q = []                                      # thương lưu sẵn: không phụ thuộc f, beta
for i in range(NV):
    qi, ri = (lagrange_basis(F, NV, i) * (T - t[i])).divmod(ZV)
    assert ri.is_zero()
    Q.append(qi)
QA_cache = Poly.zero(F)
for i in nz:
    QA_cache = QA_cache + Q[i] * A_vals[i]
print("Q_A = tổng A_i·Q_i trên %d chỉ số khác 0: %s" % (len(nz), QA_cache == QA))

srs = kzg.setup(2 * NV, tau=123)
with oc.scope():
    cmQ = [kzg.commit(srs, qi) for qi in Q]  # tiền xử lý: |V| cam kết
with oc.scope() as c1:
    cm_truc_tiep = kzg.commit(srs, QA)
with oc.scope() as c2:
    cm_cache = oc.msm([cmQ[i] for i in nz], [A_vals[i] for i in nz])
print("[Q_A]_1 trực tiếp: MSM cỡ %d; từ thương lưu sẵn: MSM cỡ %d; bằng nhau: %s"
      % (c1["msm_size"], c2["msm_size"], cm_truc_tiep == cm_cache))

# (4) Caulk: thương của tập con từ các bằng chứng mở tính sẵn
I = nz
ZI = Poly.const(F, 1)
for i in I:
    ZI = ZI * Poly(F, [-V[i], 1])
CI = lagrange_interpolate(F, [V[i] for i in I], [t[i] for i in I])
HI, du2 = (T - CI).divmod(ZI)
W = {i: (T - t[i]).divmod(Poly(F, [-V[i], 1]))[0] for i in I}
mu = {}
for i in I:
    d = F(1)
    for j in I:
        if j != i:
            d = d * (V[i] - V[j])
    mu[i] = F(1) / d
HI_cache = Poly.zero(F)
for i in I:
    HI_cache = HI_cache + W[i] * mu[i]
print("Caulk: I = %s, T - C_I chia hết cho Z_I: %s; bậc H_I = %d; H_I = tổng mu_i·W_i: %s"
      % (I, du2.is_zero(), HI.degree(), HI_cache == HI))
print("  mu_i =", [int(mu[i]) for i in I], "; C_I(omega^i) =", [int(CI(V[i])) for i in I])

# (5) Lasso: phân rã bảng 10 bit thành hai bảng con 5 bit
def bits5(x):
    return [(x >> (4 - s)) & 1 for s in range(5)]        # x_1 là bit cao
mle_id = lambda xs: sum((F(2 ** (4 - s)) * F(xs[s]) for s in range(5)), F(0))
for v in (17, 23, 391, 840 % 1024):
    r1, r2 = v >> 5, v & 31
    print("Lasso: %4d = 32·%2d + %2d -> T1[r1] = %2d, T2[r2] = %2d (MLE trên bit: %d, %d)"
          % (v, r1, r2, r1, r2, int(mle_id(bits5(r1))), int(mle_id(bits5(r2)))))
    assert 32 * r1 + r2 == v and int(mle_id(bits5(r1))) == r1
n_tra, c, d = 2, 2, 1024
print("Lasso cam kết (c = %d, alpha = %d): 3cn + alpha·d^(1/c) = %d; c = 1: n + d = %d"
      % (c, c, 3 * c * n_tra + c * 32, n_tra + d))

assert sum(A_vals, F(0)) == sum(B_vals, F(0)) and du.is_zero()
assert du2.is_zero() and HI_cache == HI
assert QA_cache == QA and cm_truc_tiep == cm_cache
