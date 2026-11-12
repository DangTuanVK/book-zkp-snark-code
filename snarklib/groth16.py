"""Groth16 (Groth, EUROCRYPT 2016, Mục 3.2) trên đường cong toy của sách.

LƯU Ý: đường cong toy y^2 = x^3 + x trên F_5387 cho pairing ĐỐI XỨNG (Type 1):
G1 = G2 = <G> và e(P, Q) = e(Q, P). Để mã đọc giống văn bản viết theo Type 3,
API vẫn tách hai phần tử sinh g1 (cho [·]_1) và g2 (cho [·]_2); trên toy cả hai
đều là G1 = (2658, 1349). Kích thước r = 449 chỉ để tính tay, không an toàn.

Quy trình:
    qap = r1cs_to_qap(r1cs)                 # nội suy trên H = <ω_n> ⊂ F_449^*
    pk, vk, tau = setup(qap, g1, g2, rng)   # tau = (α, β, γ, δ, x): "chất thải độc hại"
    proof = prove(pk, qap, z, rng)          # z = (1, x, w)
    verify(vk, x, proof)                    # e(A,B) = e(α,β)·e(IC,γ)·e(C,δ)

Ký hiệu theo bài báo: a_0 = 1, a_1..a_ℓ công khai, a_{ℓ+1}..a_m bí mật;
u_i, v_i, w_i là đa thức cột i của A, B, C; t(X) = X^n - 1.
"""

import random
from collections import namedtuple

from .field import GF
from .poly import Poly, vanishing_poly
from .fft import intt
from .pairing import e
from .curve import R_TOY

QAP = namedtuple("QAP", "F n ell u v w t m_r1cs")
ProvingKey = namedtuple("ProvingKey", "alpha1 beta1 delta1 beta2 delta2 xs1 xs2 K1 H1 ell")
VerifyingKey = namedtuple("VerifyingKey", "alpha1 beta2 gamma2 delta2 IC1 e_alpha_beta")
Proof = namedtuple("Proof", "A B C")      # A, C ∈ G1; B ∈ G2


def domain_size(m):
    n = 1
    while n < m:
        n *= 2
    return n


def r1cs_to_qap(r1cs, input_consistency=True):
    """R1CS -> QAP bằng nội suy cột trên miền con nhân H, |H| = n = 2^j.

    input_consistency=True thêm ℓ+1 hàng "a_i · 0 = 0" (A có 1 ở cột i, i = 0..ℓ)
    như ark-groth16 (LibsnarkReduction): bảo đảm các u_i của phần công khai độc
    lập tuyến tính. Hàng thêm vào đứng SAU m hàng của R1CS; các hàng còn lại tới
    n là hàng 0 (giá trị 0 trên H).
    """
    F = r1cs.F
    A, B, C = [list(map(list, M)) for M in (r1cs.A, r1cs.B, r1cs.C)]
    nv = r1cs.n
    ell = r1cs.n_public
    if input_consistency:
        for i in range(ell + 1):
            A.append([1 if j == i else 0 for j in range(nv)])
            B.append([0] * nv)
            C.append([0] * nv)
    n = domain_size(len(A))

    def cols(M):
        out = []
        for i in range(nv):
            vals = [M[j][i] if j < len(M) else 0 for j in range(n)]
            out.append(Poly(F, intt(F, vals)))   # giá trị trên H -> hệ số
        return out

    return QAP(F, n, ell, cols(A), cols(B), cols(C), vanishing_poly(F, n), r1cs.m)


def combine(qap, z):
    """A(X) = Σ a_i u_i(X), B(X), C(X) và h(X) = (A·B - C)/t."""
    F = qap.F
    Ap, Bp, Cp = Poly.zero(F), Poly.zero(F), Poly.zero(F)
    for ai, ui, vi, wi in zip(z, qap.u, qap.v, qap.w):
        ai = F(ai)
        if ai != 0:
            Ap, Bp, Cp = Ap + ui * ai, Bp + vi * ai, Cp + wi * ai
    h, rem = (Ap * Bp - Cp).divmod(qap.t)
    return Ap, Bp, Cp, h, rem


def _rand_nonzero(F, rng):
    return F(rng.randrange(1, F.p))


def setup(qap, g1, g2, rng=None, tau=None):
    """Sinh (pk, vk) theo Groth16 Mục 3.2; trả thêm tau để minh họa/kiểm thử."""
    F = qap.F
    rng = rng or random.Random()
    if tau is None:
        while True:
            al, be, ga, de, x = (_rand_nonzero(F, rng) for _ in range(5))
            if qap.t(x) != 0:          # x ∉ H, nếu không t(x) = 0
                break
    else:
        al, be, ga, de, x = (F(v) for v in tau)
    n, ell = qap.n, qap.ell
    tx = qap.t(x)
    gi, di = ga.inverse(), de.inverse()
    lin = [be * u(x) + al * v(x) + w(x) for u, v, w in zip(qap.u, qap.v, qap.w)]
    xp = [F(1)]
    for _ in range(n - 1):
        xp.append(xp[-1] * x)
    pk = ProvingKey(
        alpha1=int(al) * g1, beta1=int(be) * g1, delta1=int(de) * g1,
        beta2=int(be) * g2, delta2=int(de) * g2,
        xs1=[int(t) * g1 for t in xp], xs2=[int(t) * g2 for t in xp],
        K1=[int(l * di) * g1 for l in lin[ell + 1:]],          # i = ℓ+1..m
        H1=[int(xp[i] * tx * di) * g1 for i in range(n - 1)],  # i = 0..n-2
        ell=ell)
    alpha1 = int(al) * g1
    beta2 = int(be) * g2
    vk = VerifyingKey(alpha1=alpha1, beta2=beta2, gamma2=int(ga) * g2,
                      delta2=int(de) * g2,
                      IC1=[int(l * gi) * g1 for l in lin[:ell + 1]],   # i = 0..ℓ
                      e_alpha_beta=e(alpha1, beta2))
    return pk, vk, (al, be, ga, de, x)


def msm(scalars, points):
    """Σ s_i · P_i (nhân vô hướng đơn giản, đủ cho kích thước toy)."""
    acc = points[0].E.infinity()
    for s, P in zip(scalars, points):
        s = int(s)
        if s:
            acc = acc + s * P
    return acc


def prove(pk, qap, z, rng=None, r=None, s=None):
    """π = ([A]_1, [B]_2, [C]_1) với r, s ngẫu nhiên (cho trước được để tái lập)."""
    F = qap.F
    rng = rng or random.Random()
    r = F(rng.randrange(F.p)) if r is None else F(r)
    s = F(rng.randrange(F.p)) if s is None else F(s)
    Ap, Bp, _, h, rem = combine(qap, z)
    if not rem.is_zero():
        raise ValueError("witness không thỏa R1CS: t(X) không chia hết A·B - C")
    ell = pk.ell
    A1 = pk.alpha1 + msm(Ap.coeffs, pk.xs1) + int(r) * pk.delta1
    B2 = pk.beta2 + msm(Bp.coeffs, pk.xs2) + int(s) * pk.delta2
    B1 = pk.beta1 + msm(Bp.coeffs, pk.xs1) + int(s) * pk.delta1
    C1 = (msm(z[ell + 1:], pk.K1) + msm(h.coeffs, pk.H1)
          + int(s) * A1 + int(r) * B1 - int(r * s) * pk.delta1)
    return Proof(A1, B2, C1)


def public_commitment(vk, x):
    """IC = [Σ_{i=0..ℓ} a_i (βu_i + αv_i + w_i)(x)/γ]_1 với a_0 = 1."""
    return msm([1] + [int(v) for v in x], vk.IC1)


def verify(vk, x, proof):
    """Kiểm e(A, B) = e(α, β) · e(IC, γ) · e(C, δ)."""
    if len(x) != len(vk.IC1) - 1:
        raise ValueError("số đầu vào công khai không khớp")
    lhs = e(proof.A, proof.B)
    rhs = vk.e_alpha_beta * e(public_commitment(vk, x), vk.gamma2) * e(proof.C, vk.delta2)
    return lhs == rhs


def rerandomize(vk, proof, r, s):
    """Tái ngẫu nhiên hóa: A' = r^{-1}A, B' = rB + rs·[δ]_2, C' = C + s·A.

    Kiểm: A'·B' = AB + sδA = (C + sA)δ - Cδ + AB, nên phương trình verify giữ nguyên.
    Không cần witness, chỉ cần vk (để có [δ]_2).
    """
    F = GF(R_TOY)
    r, s = F(r), F(s)
    if r == 0:
        raise ValueError("r phải khác 0")
    return Proof(int(r.inverse()) * proof.A,
                 int(r) * proof.B + int(r * s) * vk.delta2,
                 proof.C + int(s) * proof.A)


def simulate(qap, tau, x, g1, g2, rng=None):
    """Sim(τ, x): A, B ← F ngẫu nhiên; C = (AB - αβ - Σ_{i≤ℓ} a_i(βu_i+αv_i+w_i)(x))/δ."""
    F = qap.F
    rng = rng or random.Random()
    al, be, ga, de, xx = (F(v) for v in tau)
    A, B = F(rng.randrange(F.p)), F(rng.randrange(F.p))
    a = [F(1)] + [F(v) for v in x]
    pub = sum((ai * (be * qap.u[i](xx) + al * qap.v[i](xx) + qap.w[i](xx))
               for i, ai in enumerate(a)), F(0))
    C = (A * B - al * be - pub) / de
    return Proof(int(A) * g1, int(B) * g2, int(C) * g1)
