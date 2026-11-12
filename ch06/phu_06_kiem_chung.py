"""Chương 6 -- Kiểm chứng phụ cho Mục 6.1-6.2 (không in trong sách, chỉ sinh số liệu).

Các con số dùng ở phần bảng cổng, vòng 1-3 và bổ đề hoán vị:
  (1) ba coset H, 2H, 3H phân biệt; nhãn σ*(a0) = k2·ω^19;
  (2) làm mù: a(ω^19) = 17 dù a(X) có thêm (b1 X + b2) Z_H(X);
  (3) ba tỉ số đầu của z(X) với β = 216, γ = 28 (seed 2019 như Mã nguồn 6.3);
  (4) bậc của ba đa thức ràng buộc và kiểm tra tại ζ: tử số(ζ) = t(ζ)·Z_H(ζ);
  (5) vét cạn mọi (β, γ) ∈ F_449^2: witness trung thực và witness sai dây
      (a19 = 18, c19 = 414) cho tích F(β, γ) = G(β, γ) bao nhiêu lần;
  (6) xác suất Prover phải dừng vì mẫu số của z bằng 0.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib import plonk
from snarklib.plonk import F, pi_poly
from snarklib.poly import Poly
from snarklib.dsl import mach_phan_tich

c = mach_phan_tich()[0]
circ = c.compile_plonkish()
a, b, cc, x = c.plonkish_assignment({"N": 391, "a": 17, "b": 23})
pk, vk = plonk.setup(circ, tau=123)
n, om, k1, k2, Q = pk.n, pk.omega, pk.k1, pk.k2, pk.polys
p = 449

# (1) coset và nhãn
print("(1) 2^32 = %d, 3^32 = %d, (2/3)^32 = %d (đều khác 1); số coset của H trong F*: %d"
      % (int(F(2) ** n), int(F(3) ** n), int((F(2) / F(3)) ** n), (p - 1) // n))
print("    ω^j, j = 0..31:", [int(om ** j) for j in range(n)])
print("    k2·ω^19 = %d = σ*(a0);  Z_H(τ = 123) = %d" % (int(k2 * om ** 19), int(F(123) ** n - 1)))

P = plonk.Prover(pk, a, b, cc, x, random.Random(2019))
P.run()
be, ga, al, ze = P.beta, P.gamma, P.alpha, P.zeta

# (2) làm mù
print("(2) b1..b11 =", [int(t) for t in P.b[1:]])
print("    a(ω^19) = %d, a(τ) = %d; bậc a, z = %d, %d"
      % (int(P.a_p(om ** 19)), int(P.a_p(F(123))), P.a_p.degree(), P.z_p.degree()))

# (3) ba tỉ số đầu của z
w, ids, ss = P.w, pk.ids, pk.sstar
for i in range(3):
    num = [w[j][i] + be * ids[j * n + i] + ga for j in range(3)]
    den = [w[j][i] + be * ss[j * n + i] + ga for j in range(3)]
    r = num[0] * num[1] * num[2] / (den[0] * den[1] * den[2])
    print("(3) hàng %d: nhãn %s, σ* %s, tử %s, mẫu %s, tỉ số %d, z(ω^%d) = %d"
          % (i, [int(ids[j * n + i]) for j in range(3)], [int(ss[j * n + i]) for j in range(3)],
             [int(t) for t in num], [int(t) for t in den], int(r), i + 1, int(P.z_vals[i + 1])))

# (4) ba ràng buộc và kiểm tra tại ζ
X = Poly.X(F)
A, B, C, Z = P.a_p, P.b_p, P.c_p, P.z_p
gate = A * B * Q["qM"] + A * Q["qL"] + B * Q["qR"] + C * Q["qO"] + pi_poly(n, x) + Q["qC"]
perm = ((A + X * be + ga) * (B + X * (be * k1) + ga) * (C + X * (be * k2) + ga) * Z
        - (A + Q["s1"] * be + ga) * (B + Q["s2"] * be + ga) * (C + Q["s3"] * be + ga) * Z.scale(om))
init = (Z - 1) * Q["L0"]
for ten, f in (("cổng", gate), ("hoán vị", perm), ("khởi đầu", init)):
    print("(4) ràng buộc %-8s bậc %3d, chia hết cho Z_H: %s, giá trị tại ζ = %d"
          % (ten, f.degree(), f.divmod(P.ZH)[1].is_zero(), int(f(ze))))
tong = gate(ze) + al * perm(ze) + al * al * init(ze)
print("    α = %d, ζ = %d: tử số(ζ) = %d; t(ζ) = %d, Z_H(ζ) = %d, t(ζ)·Z_H(ζ) = %d; bậc t = %d"
      % (int(al), int(ze), int(tong), int(P.t_p(ze)), int(ze ** n - 1),
         int(P.t_p(ze) * (ze ** n - 1)), P.t_p.degree()))
print("    t_lo(ζ) + ζ^n t_mid(ζ) + ζ^2n t_hi(ζ) = %d"
      % int(P.t_lo(ze) + ze ** n * P.t_mid(ze) + ze ** (2 * n) * P.t_hi(ze)))


# (5) vét cạn (β, γ): chỉ các ô không cố định của σ mới đóng góp vào F/G
def vet_can(a, b, c):
    ww = [[int(t) for t in v] + [0] * (n - len(v)) for v in (a, b, c)]
    I, S = [int(v) for v in ids], [int(v) for v in ss]
    cells = [(ww[j][i], I[j * n + i], S[j * n + i]) for i in range(n) for j in range(3)
             if I[j * n + i] != S[j * n + i]]
    bang = khong = 0
    for Be in range(p):
        for Ga in range(p):
            nu = de = 1
            for v, u, s in cells:
                nu = nu * (v + Be * u + Ga) % p
                de = de * (v + Be * s + Ga) % p
            if de == 0:
                khong += 1
            elif nu == de:
                bang += 1
    return len(cells), bang, khong


a2, c2 = list(a), list(cc)
a2[19], c2[19] = 18, 18 * 23 % p
for ten, ww in (("trung thực", (a, b, cc)), ("sai dây a19 = 18, c19 = 414", (a2, b, c2))):
    m, bang, khong = vet_can(*ww)
    print("(5) %-28s: %d ô không cố định; trong %d cặp (β, γ): mẫu = 0: %d, F = G (mẫu ≠ 0): %d"
          % (ten, m, p * p, khong, bang))
print("    cận Schwartz-Zippel 58/449 = %.4f" % (58 / p))

# (6) xác suất dừng
print("(6) 1 - (1 - 1/449)^93 = %.4f; 3n/r = %.4f" % (1 - (1 - 1 / p) ** (3 * (n - 1)), 3 * n / p))
assert int(P.z_vals[1]) == 151 and int(tong) == int(P.t_p(ze) * (ze ** n - 1))

# (7) cam kết [a]_1 = a(τ)·G1 và các số hạng sai số ở kích thước thật (BN254, n = 2^20)
import math
print("(7) a(τ)·G1 = %s, [a]_1 trong proof = %s"
      % (int(P.a_p(F(123))) * pk.srs.g1[0], P.cm_abc[0]))
r_bn = 21888242871839275222246405745257275088548364400416034343698204186575808495617
for ten, d in (("3(n-1) (mẫu số bằng 0)", 3 * (n - 1)), ("2 (gộp bằng α)", 2),
               ("3n (Schwartz-Zippel, hoán vị)", 3 * n), ("n (ζ ∈ H)", n)):
    print("    %-30s toy n = 32: %.4f | BN254, n = 2^20: 2^%.1f"
          % (ten, d / p, math.log2((d if d < 32 else d * 2 ** 20 // 32) / r_bn)))
