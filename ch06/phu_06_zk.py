"""Chương 6 -- Kiểm chứng phụ cho phần zero-knowledge (prop:ch6-zk), không in trong sách.

Mạch N = 391 có đúng 4 witness: (17, 23), (23, 17), (28, 30), (30, 28).
  (1) Không làm mù: các giá trị a(ζ), b(ζ), c(ζ), z(ζω) tất định theo witness,
      bốn witness cho bốn bộ khác nhau -> giá trị mở để lộ witness.
  (2) Có làm mù: với cửa sập τ (vai trò bộ mô phỏng), giải ρ'_1..ρ'_11 để Prover
      dùng witness KHÁC tạo ra proof GIỐNG HỆT từng byte proof của (17, 23).
      Hệ phương trình: 2 ẩn cho mỗi a, b, c (điểm τ, ζ), 3 ẩn cho z (điểm τ, ωτ, ζω),
      2 ẩn cho t_lo, t_mid (điểm τ).
  (3) Nếu z chỉ có 2 hệ số làm mù: hệ 3 phương trình 2 ẩn vô nghiệm -> phân biệt được.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib import plonk
from snarklib.plonk import F
from snarklib.poly import Poly
from snarklib.fft import intt
from snarklib.dsl import mach_phan_tich

c = mach_phan_tich()[0]
circ = c.compile_plonkish()
pk, vk = plonk.setup(circ, tau=123)
n, om, tau = pk.n, pk.omega, F(123)
ZH = lambda X: X ** n - 1
witness = [(17, 23), (23, 17), (28, 30), (30, 28)]
cols = {w: c.plonkish_assignment({"N": 391, "a": w[0], "b": w[1]}) for w in witness}

P0 = plonk.Prover(pk, *cols[(17, 23)], random.Random(2019))
proof0 = P0.run()
be, ga, ze = P0.beta, P0.gamma, P0.zeta


def interp(vals):
    v = [F(t) for t in vals] + [F(0)] * (n - len(vals))
    return Poly(F, intt(F, v))


def z_values(a, b, cc):
    w = [[F(t) for t in col] + [F(0)] * (n - len(col)) for col in (a, b, cc)]
    zv = [F(1)]
    for i in range(n - 1):
        num = den = F(1)
        for j in range(3):
            num *= w[j][i] + be * pk.ids[j * n + i] + ga
            den *= w[j][i] + be * pk.sstar[j * n + i] + ga
        zv.append(zv[-1] * num / den)
    return zv


# (1) không làm mù
print("(1) ζ = %d; giá trị mở nếu KHÔNG làm mù (β, γ, ζ lấy từ transcript seed 2019):" % int(ze))
bo = set()
for w in witness:
    a, b, cc, _ = cols[w]
    v = tuple(int(interp(col)(ze)) for col in (a, b, cc)) + (int(interp(z_values(a, b, cc))(ze * om)),)
    bo.add(v)
    print("    witness %-9s a(ζ), b(ζ), c(ζ), z(ζω) = %s" % (w, v))
print("    số bộ khác nhau: %d / %d" % (len(bo), len(witness)))


def solve(rows, rhs):
    """Giải hệ tuyến tính nhỏ trên F bằng khử Gauss; trả None nếu vô nghiệm."""
    m = [list(r) + [y] for r, y in zip(rows, rhs)]
    k, piv = len(rows[0]), []
    rr = 0
    for col in range(k):
        p = next((i for i in range(rr, len(m)) if m[i][col] != 0), None)
        if p is None:
            continue
        m[rr], m[p] = m[p], m[rr]
        inv = m[rr][col].inverse()
        m[rr] = [t * inv for t in m[rr]]
        for i in range(len(m)):
            if i != rr and m[i][col] != 0:
                f = m[i][col]
                m[i] = [s - f * t for s, t in zip(m[i], m[rr])]
        piv.append(col)
        rr += 1
    if any(all(t == 0 for t in row[:-1]) and row[-1] != 0 for row in m):
        return None
    x = [F(0)] * k
    for i, col in enumerate(piv):
        x[col] = m[i][-1]
    return x


class FixedRng:
    def __init__(self, vals):
        self.vals = [int(v) for v in vals]

    def randrange(self, _):
        return self.vals.pop(0)


# (2) giả lập: witness khác, cùng proof
for w in witness[1:]:
    a, b, cc, x = cols[w]
    rho = []
    for k, col in enumerate((a, b, cc)):
        f0 = (P0.a_p, P0.b_p, P0.c_p)[k]
        g = interp(col)
        pts = [tau, ze]
        rho += solve([[X * ZH(X), ZH(X)] for X in pts], [f0(X) - g(X) for X in pts])
    zl = interp(z_values(a, b, cc))
    pts = [tau, om * tau, ze * om]
    rho += solve([[X * X * ZH(X), X * ZH(X), ZH(X)] for X in pts],
                 [P0.z_p(X) - zl(X) for X in pts])
    thu = plonk.Prover(pk, a, b, cc, x, FixedRng(rho + [0, 0]))
    thu.run()
    r10 = (P0.t_lo(tau) - thu.t_lo(tau)) / tau ** n
    r11 = (P0.t_mid(tau) - thu.t_mid(tau) + r10) / tau ** n
    P1 = plonk.Prover(pk, a, b, cc, x, FixedRng(rho + [r10, r11]))
    proof1 = P1.run()
    print("(2) witness %-9s ρ' = %s -> proof trùng proof của (17, 23): %s, verify: %s"
          % (w, [int(t) for t in rho + [r10, r11]], proof1 == proof0, plonk.verify(vk, x, proof1)))
    assert proof1 == proof0

# (3) z chỉ có 2 hệ số làm mù
a, b, cc, _ = cols[(23, 17)]
zl = interp(z_values(a, b, cc))
pts = [tau, om * tau, ze * om]
kq = solve([[X * ZH(X), ZH(X)] for X in pts], [P0.z_p(X) - zl(X) for X in pts])
print("(3) z với 2 hệ số làm mù: hệ 3 phương trình (τ, ωτ, ζω) 2 ẩn có nghiệm? %s" % (kq is not None))
print("    số điểm lộ của mỗi đa thức: a, b, c: 2 (τ, ζ); z: 3 (τ, ωτ qua t(τ), ζω); "
      "t_lo, t_mid, t_hi: 3 cam kết, 1 ràng buộc tổng")
