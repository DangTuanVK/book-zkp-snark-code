"""Chương 5 (mã phụ) -- QAP của mạch N = 391 nhìn như một PCP tuyến tính.

Số liệu cho Mục 5.2 (linear PCP, LIP, NILP):
  (1) Proof tuyến tính pi = (z, hệ số của h); bốn truy vấn q1..q4 tại điểm r;
      Verifier bậc 2: <q1,z>·<q2,z> - <q3,z> = <q4,h>·t(r).
  (2) Cận Schwartz--Zippel trên trường toy: Prover dùng MỘT z sai (N = 389)
      và thương h của phép chia có dư; đếm số r làm phép kiểm lọt.
  (3) Prover trả lời mỗi truy vấn bằng một hàm tuyến tính khác nhau (zA, zB, zC)
      chứng minh được N = 389 (không có witness); truy vấn kiểm nhất quán
      q5 = b1·q1 + b2·q2 + b3·q3 (LIP) bắt được.
Quy ước z = (1, N, a, b, a0..a4, b0..b4); QAP có 13 + 2 hàng, |H| = 16.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.field import GF
from snarklib.dsl import mach_phan_tich
from snarklib.groth16 import r1cs_to_qap, combine

F = GF(449)
c = mach_phan_tich()[0]
r1cs = c.compile_r1cs()
qap = r1cs_to_qap(r1cs)
n = qap.n
x, w = c.r1cs_assignment({"N": 391, "a": 17, "b": 23})
z = [F(1)] + list(x) + list(w)
A, B, C, h, rem = combine(qap, z)
assert rem.is_zero() if hasattr(rem, "is_zero") else all(v == 0 for v in rem.coeffs)

r = F(427)                           # trùng với x của CRS Groth16 (seed 2016) ở Mục 5.3
q1 = [u(r) for u in qap.u]
q2 = [v(r) for v in qap.v]
q3 = [wp(r) for wp in qap.w]
hc = list(h.coeffs) + [F(0)] * (n - 1 - len(h.coeffs))
q4 = [r ** k for k in range(n - 1)]
dot = lambda q, v: sum((qi * vi for qi, vi in zip(q, v)), F(0))
a1, a2, a3, a4 = dot(q1, z), dot(q2, z), dot(q3, z), dot(q4, hc)
tr = qap.t(r)
print("== (1) QAP như PCP tuyến tính, N = 391, (a, b) = (17, 23) ==")
print("|z| = %d, deg h = %d -> proof pi = (z, h) dài %d phần tử; |H| = %d, t(X) = X^%d - 1"
      % (len(z), h.degree(), len(z) + len(hc), n, n))
print("r = %d: <q1,z> = A(r) = %d, <q2,z> = B(r) = %d, <q3,z> = C(r) = %d, <q4,h> = h(r) = %d, t(r) = %d"
      % (int(r), int(a1), int(a2), int(a3), int(a4), int(tr)))
print("kiểm bậc 2: A·B - C = %d, h·t = %d -> %s" % (int(a1 * a2 - a3), int(a4 * tr), a1 * a2 - a3 == a4 * tr))
assert a1 == A(r) and a2 == B(r) and a3 == C(r) and a4 == h(r) and a1 * a2 - a3 == a4 * tr

# ---------------------------------------------------------------- (2) Schwartz--Zippel
print("== (2) Một z sai, N = 389 (thay N, giữ a, b) ==")
zb = list(z)
zb[1] = F(389)
Ab, Bb, Cb, hb, remb = combine(qap, zb)
deg_rem = remb.degree()
pass_r = [s for s in range(449) if remb(F(s)) == 0]
H = {int(F(349) ** k) for k in range(n)}          # H = <omega_16>, omega_16 = 349
ngoai_H = [s for s in pass_r if s not in H]
thieu = sorted(H - set(pass_r))
print("phần dư R(X) = A·B - C - h·t có bậc %d; số r trong F_449 làm phép kiểm lọt: %d %s"
      % (deg_rem, len(pass_r), pass_r))
print("trong đó thuộc H: %d (mọi điểm của H trừ %s = omega^12, hàng cổng nhân); ngoài H: %d"
      % (len(pass_r) - len(ngoai_H), thieu, len(ngoai_H)))
assert thieu == [int(F(349) ** 12)] and not ngoai_H
print("cận tổng quát: deg(A·B - C - h·t) <= 2(n-1) = %d -> r đều trên F \\ H: xác suất lọt <= %d/%d = %.4f"
      % (2 * (n - 1), 2 * (n - 1), 449 - n, 2 * (n - 1) / (449 - n)))
wit389 = [(a, b) for a in range(32) for b in range(32) if (a * b) % 449 == 389]
print("vét cạn a, b trong [0, 31] với a·b = 389 (mod 449):", wit389)

# ---------------------------------------------------------------- (3) LIP
print("== (3) Ba hàm tuyến tính khác nhau cho ba truy vấn, N = 389 ==")
idx = {nm: j for j, nm in enumerate(["1", "N", "a", "b"] + ["a%d" % i for i in range(5)]
                                      + ["b%d" % i for i in range(5)])}
zC = list(zb)                       # a = 17, b = 23 cùng các bit đúng, N = 389
zA = list(zC); zA[idx["a"]] = F(1)  # chỉ đổi a ở hàm trả lời q1
zB = list(zC); zB[idx["b"]] = F(389)  # chỉ đổi b ở hàm trả lời q2
AA = sum((qap.u[i] * zA[i] for i in range(len(z)) if zA[i] != 0), type(A).zero(F))
BB = sum((qap.v[i] * zB[i] for i in range(len(z)) if zB[i] != 0), type(A).zero(F))
CC = sum((qap.w[i] * zC[i] for i in range(len(z)) if zC[i] != 0), type(A).zero(F))
hh, rr = (AA * BB - CC).divmod(qap.t)
ok_poly = all(v == 0 for v in rr.coeffs) if hasattr(rr, "coeffs") else rr.is_zero()
b1, b2, b3 = dot(q1, zA), dot(q2, zB), dot(q3, zC)
print("zA: a = 1; zB: b = 389; zC: a = 17, b = 23; mọi zX có N = 389")
print("A_zA·B_zB - C_zC chia hết cho t(X): %s; tại r: %d·%d - %d = %d, h(r)·t(r) = %d"
      % (ok_poly, int(b1), int(b2), int(b3), int(b1 * b2 - b3), int(hh(r) * tr)))
assert ok_poly and b1 * b2 - b3 == hh(r) * tr
be = [F(2), F(3), F(5)]             # hệ số bí mật của truy vấn kiểm nhất quán
q5 = [be[0] * u + be[1] * v + be[2] * wq for u, v, wq in zip(q1, q2, q3)]
for nm, zz in (("trung thực (một z)", z), ("gian lận (zA, zB, zC)", None)):
    if zz is not None:
        ans = (dot(q1, zz), dot(q2, zz), dot(q3, zz))
        a5 = dot(q5, zz)
    else:
        ans = (b1, b2, b3)
        a5 = dot(q5, zC)            # q5 trả lời bằng một trong ba hàm (ở đây zC)
    ok5 = a5 == be[0] * ans[0] + be[1] * ans[1] + be[2] * ans[2]
    print("q5 = 2·q1 + 3·q2 + 5·q3, Prover %s: <q5,.> = %d, 2·a1 + 3·a2 + 5·a3 = %d -> %s"
          % (nm, int(a5), int(be[0] * ans[0] + be[1] * ans[1] + be[2] * ans[2]), ok5))
