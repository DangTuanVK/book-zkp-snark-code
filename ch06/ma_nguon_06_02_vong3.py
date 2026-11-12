"""Chương 6 -- Vòng 3: tính t(X) bằng coset FFT (lst:ch6-round3).

F_449 chỉ có căn đơn vị bậc ≤ 64, không có miền cỡ 4n = 128. Ta đánh giá tử
số trên 4 coset s_j·H (mỗi coset một NTT cỡ n = 32), chia điểm-điểm cho
Z_H(s_j ω^i) = s_j^n - 1 (hằng trên coset), sau đó ghép bằng phần dư Trung Hoa.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib import plonk, opcount as oc
from snarklib.plonk import F, eval_on_coset, interpolate_from_cosets, pi_poly
from snarklib.dsl import mach_phan_tich

c = mach_phan_tich()[0]
circ = c.compile_plonkish()
a, b, cc, x = c.plonkish_assignment({"N": 391, "a": 17, "b": 23})
pk, vk = plonk.setup(circ, tau=123)
P = plonk.Prover(pk, a, b, cc, x, random.Random(6))
P.run()                                     # lấy a(X), b(X), c(X), z(X), β, γ, α thật
n, om, be, ga, al = pk.n, pk.omega, P.beta, P.gamma, P.alpha
k1, k2, Q = pk.k1, pk.k2, pk.polys

def tu_so(s):
    """Tử số của t(X) tại các điểm s·ω^i của coset s·H."""
    A, B, C, Z = (eval_on_coset(f, s, n) for f in (P.a_p, P.b_p, P.c_p, P.z_p))
    Zw = Z[1:] + Z[:1]                      # z(ω·sω^i) = z(sω^{i+1}): xoay 1 vị trí
    q = {k: eval_on_coset(Q[k], s, n) for k in ("qM", "qL", "qR", "qO", "qC",
                                              "s1", "s2", "s3", "L0")}
    PI = eval_on_coset(pi_poly(n, x), s, n)
    out = []
    for i in range(n):
        X = s * om ** i
        cong = A[i]*B[i]*q["qM"][i] + A[i]*q["qL"][i] + B[i]*q["qR"][i] \
            + C[i]*q["qO"][i] + PI[i] + q["qC"][i]
        hv1 = (A[i] + be*X + ga) * (B[i] + be*k1*X + ga) * (C[i] + be*k2*X + ga) * Z[i]
        hv2 = (A[i] + be*q["s1"][i] + ga) * (B[i] + be*q["s2"][i] + ga) \
            * (C[i] + be*q["s3"][i] + ga) * Zw[i]
        out.append(cong + al*(hv1 - hv2) + al*al*(Z[i] - 1)*q["L0"][i])
    return out

with oc.scope() as dem:
    t_vals = []
    for s in pk.shifts:
        zh_inv = (s ** n - 1).inverse()     # Z_H hằng trên coset s·H
        t_vals.append([v * zh_inv for v in tu_so(s)])
    t = interpolate_from_cosets(t_vals, pk.shifts, n)

t_chia, du = P.quotient_naive()             # đối chiếu: chia đa thức trực tiếp
print("coset s_j =", [int(s) for s in pk.shifts], " s_j^n =", [int(s ** n) for s in pk.shifts])
print("β, γ, α =", int(be), int(ga), int(al))
print("bậc t(X) = %d (≤ 3n + 5 = %d); số hệ số khác 0 = %d" % (t.degree(), 3 * n + 5,
      sum(1 for v in t.coeffs if v != 0)))
print("t(X) 4 hệ số đầu:", [int(v) for v in t.coeffs[:4]], "... hệ số bậc 101:", int(t.coeffs[101]))
print("số NTT cỡ %d trong khối trên: %d (kể cả NTT của selector, trong Prover thật đã tiền xử lý)"
      % (n, dem["ntt"]))
print("trùng với phép chia đa thức:", t == t_chia, "| phần dư bằng 0:", du.is_zero())
assert t == t_chia and du.is_zero() and t.degree() <= 3 * n + 5
