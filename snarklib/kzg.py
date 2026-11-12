"""Cam kết đa thức KZG và các biến thể mở nhiều đa thức/nhiều điểm (Chương 3).

Ký hiệu: [a]_1 = a·G1, [a]_2 = a·G2, τ là cửa sập (toxic waste) của SRS.
Trên đường cong toy của sách, pairing là ĐỐI XỨNG (Type 1): G2 = G1 và
e = pairing Tate rút gọn trên E: y^2 = x^3 + x / F_5387, r = 449. Mọi hàm vẫn
nhận/ghi riêng phần G1 và G2 (srs.g1, srs.g2) để văn bản sách viết dạng Type 3:
thay hai điểm sinh khác nhau là chạy được trên BN254/BLS12-381.

Các biến thể:
    commit / open / verify               KZG một điểm [KZG10, §3.2]
    open_batch / verify_batch            nhiều đa thức, một điểm (γ từ Transcript)
                                         [PLONK, §3.1 bước 3]
    open_multi / verify_multi            một đa thức, nhiều điểm: Z_S và r(X)
                                         [KZG10, §3.4 "Batch Opening"]
    open_gwc / verify_gwc                nhiều đa thức, nhiều điểm: mỗi điểm một
                                         nhân chứng, gộp pairing bằng hệ số ngẫu nhiên
                                         [PLONK, §3.1, trường hợp nhiều điểm]
    open_shplonk / verify_shplonk        SHPLONK [BDFG20, ePrint 2020/081]: 2 phần
                                         tử G1 cho mọi tập điểm
    hiding_setup / commit_hiding / ...   hiding KZG kiểu PolyCommit_Ped [KZG10, §3.3]

Mọi phép toán nhóm/pairing đi qua snarklib.opcount để đếm chi phí thật.
"""

import hashlib
from collections import namedtuple

from .field import GF
from .poly import Poly, lagrange_interpolate
from .curve import toy_generator, R_TOY, COFACTOR_TOY
from .transcript import Transcript
from . import opcount as oc

R = R_TOY
F = GF(R)


class SRS:
    """SRS của KZG: g1 = [[τ^i]_1, i = 0..d], g2 = [[τ^j]_2, j = 0..d2].

    h1 (tùy chọn) = [[τ^i]H] cho hiding KZG, H có logarit rời rạc không biết.
    """

    def __init__(self, g1, g2, h1=None):
        self.g1, self.g2, self.h1 = g1, g2, h1

    @property
    def max_degree(self):
        return len(self.g1) - 1

    @property
    def G1(self):
        return self.g1[0]

    @property
    def G2(self):
        return self.g2[0]


def setup(degree, tau, g2_degree=1, G1=None, G2=None, H=None):
    """Sinh SRS bậc `degree` từ cửa sập tau (chỉ dùng cho minh họa/kiểm thử).

    Mặc định G1 là điểm sinh tất định của sách và G2 = G1 (Type 1).
    g2_degree > 1 cần cho open_multi (verifier tính [Z_S(τ)]_2).
    Thời gian sinh SRS không được đếm vào opcount.
    """
    if G1 is None:
        G1 = toy_generator()[2]
    if G2 is None:
        G2 = G1                      # toy: pairing đối xứng
    t = F(tau)
    pw = [F(1)]
    for _ in range(max(degree, g2_degree)):
        pw.append(pw[-1] * t)
    g1 = [int(pw[i]) * G1 for i in range(degree + 1)]
    g2 = [int(pw[j]) * G2 for j in range(g2_degree + 1)]
    h1 = [int(pw[i]) * H for i in range(degree + 1)] if H is not None else None
    return SRS(g1, g2, h1)


def _poly(f):
    return f if isinstance(f, Poly) else Poly(F, f)


def vanishing(points):
    """Z_S(X) = Π_{s ∈ S} (X - s)."""
    Z = Poly.const(F, 1)
    for s in points:
        Z = Z * Poly(F, [-F(s), 1])
    return Z


def _exact_div(num, den):
    q, rem = num.divmod(den)
    assert rem.is_zero(), "phép chia phải chia hết"
    return q


# ---------------------------------------------------------------------------
# KZG một điểm
# ---------------------------------------------------------------------------

def commit(srs, f, bases=None):
    """C = [f(τ)]_1 = Σ f_i·[τ^i]_1 (một MSM cỡ deg f + 1)."""
    f = _poly(f)
    bases = srs.g1 if bases is None else bases
    if f.degree() > len(bases) - 1:
        raise ValueError("bậc %d vượt SRS (%d)" % (f.degree(), len(bases) - 1))
    if f.is_zero():
        return bases[0].E.infinity()
    return oc.msm(bases[:len(f.coeffs)], f.coeffs)


def open(srs, f, z):
    """Mở f tại z: trả về (y = f(z), W = [q(τ)]_1) với q = (f - y)/(X - z)."""
    f = _poly(f)
    z = F(z)
    y = f(z)
    q = _exact_div(f - y, Poly(F, [-z, 1]))
    return y, commit(srs, q)


def verify(srs, C, z, y, W):
    """Kiểm e(C - [y]_1 + z·W, [1]_2) · e(-W, [τ]_2) = 1.

    Tương đương e(C - [y]_1, [1]_2) = e(W, [τ - z]_2) nhưng chuyển z sang G1.
    """
    lhs = oc.msm([C, srs.G1, W], [1, -F(y), F(z)])
    return oc.pairing_check([(lhs, srs.g2[0]), (-W, srs.g2[1])])


# ---------------------------------------------------------------------------
# Nhiều đa thức, một điểm — tổ hợp tuyến tính ngẫu nhiên γ
# ---------------------------------------------------------------------------

def _powers(x, n):
    out = [F(1)]
    for _ in range(n - 1):
        out.append(out[-1] * x)
    return out


def _gamma_one_point(tr, Cs, z, ys):
    tr.absorb(b"C", list(Cs))
    tr.absorb(b"z", F(z))
    tr.absorb(b"y", list(ys))
    return tr.challenge(b"gamma")


def open_batch(srs, fs, Cs, z, tr):
    """Mở k đa thức tại cùng điểm z bằng MỘT nhân chứng.

    γ lấy từ Transcript SAU khi đã hấp thụ C_i, z, y_i. Trả về (ys, W) với
    W = [h(τ)]_1, h(X) = Σ γ^(i-1) (f_i(X) - y_i)/(X - z).
    """
    fs = [_poly(f) for f in fs]
    ys = [f(F(z)) for f in fs]
    g = _powers(_gamma_one_point(tr, Cs, z, ys), len(fs))
    comb = Poly.zero(F)
    for gi, f, y in zip(g, fs, ys):
        comb = comb + (f - y) * gi
    h = _exact_div(comb, Poly(F, [-F(z), 1]))
    return ys, commit(srs, h)


def verify_batch(srs, Cs, z, ys, W, tr):
    """Verifier: F = Σ γ^(i-1) C_i, v = Σ γ^(i-1) y_i, rồi kiểm như KZG một điểm."""
    g = _powers(_gamma_one_point(tr, Cs, z, ys), len(Cs))
    v = sum((gi * y for gi, y in zip(g, ys)), F(0))
    lhs = oc.msm(list(Cs) + [srs.G1, W], g + [-v, F(z)])
    return oc.pairing_check([(lhs, srs.g2[0]), (-W, srs.g2[1])])


# ---------------------------------------------------------------------------
# Một đa thức, nhiều điểm — Z_S và nội suy r(X)
# ---------------------------------------------------------------------------

def open_multi(srs, f, S):
    """Mở f tại tập S: trả về (ys, W) với W = [(f - r)/Z_S (τ)]_1, r nội suy f trên S."""
    f = _poly(f)
    ys = [f(F(s)) for s in S]
    r = lagrange_interpolate(F, S, ys)
    psi = _exact_div(f - r, vanishing(S))
    return ys, commit(srs, psi)


def verify_multi(srs, C, S, ys, W):
    """Kiểm e(C - [r(τ)]_1, [1]_2) = e(W, [Z_S(τ)]_2); cần len(srs.g2) > |S|."""
    if len(srs.g2) <= len(S):
        raise ValueError("SRS cần [τ^j]_2 tới j = |S| = %d" % len(S))
    r = lagrange_interpolate(F, S, ys)
    ZS = vanishing(S)
    lhs = oc.sub(C, commit(srs, r))
    ZS2 = oc.msm(srs.g2[:len(ZS.coeffs)], ZS.coeffs, group="G2")
    return oc.pairing_check([(lhs, srs.g2[0]), (-W, ZS2)])


# ---------------------------------------------------------------------------
# Nhiều đa thức, nhiều điểm (cách của PLONK §3.1): mỗi điểm phân biệt một nhân
# chứng; verifier gộp các phương trình bằng hệ số ngẫu nhiên r'.
# ---------------------------------------------------------------------------

GWCProof = namedtuple("GWCProof", "ys Ws")


def _distinct(Ss):
    T = []
    for S in Ss:
        for s in S:
            if F(s) not in T:
                T.append(F(s))
    return T


def _absorb_claims(tr, Cs, Ss, ys):
    tr.absorb(b"C", list(Cs))
    tr.absorb(b"S", [[F(s) for s in S] for S in Ss])
    tr.absorb(b"y", [list(y) for y in ys])


def open_gwc(srs, fs, Cs, Ss, tr):
    """fs[i] mở tại tập Ss[i]. Trả về GWCProof(ys, Ws): một W cho mỗi điểm phân biệt."""
    fs = [_poly(f) for f in fs]
    ys = [[f(F(s)) for s in S] for f, S in zip(fs, Ss)]
    _absorb_claims(tr, Cs, Ss, ys)
    Ws = []
    for t in _distinct(Ss):
        gam = tr.challenge(b"gamma")
        comb, gi = Poly.zero(F), F(1)
        for f, S, y in zip(fs, Ss, ys):
            for s, v in zip(S, y):
                if F(s) == t:
                    comb = comb + (f - v) * gi
                    gi = gi * gam
        Ws.append(commit(srs, _exact_div(comb, Poly(F, [-t, 1]))))
    return GWCProof(ys, Ws)


def verify_gwc(srs, Cs, Ss, proof, tr):
    """Kiểm e(Σ_j r'^j (F_j - [v_j]_1 + t_j W_j), [1]_2) · e(-Σ_j r'^j W_j, [τ]_2) = 1."""
    _absorb_claims(tr, Cs, Ss, proof.ys)
    T = _distinct(Ss)
    gams = [tr.challenge(b"gamma") for _ in T]
    tr.absorb(b"W", list(proof.Ws))
    rr = _powers(tr.challenge(b"r'"), len(T))
    # gom toàn bộ vế trái vào một MSM trên (C_1..C_k, G1, W_1..W_t*)
    sc = [F(0)] * len(Cs)
    sG = F(0)
    for j, t in enumerate(T):
        gi = F(1)
        for i, (S, y) in enumerate(zip(Ss, proof.ys)):
            for s, v in zip(S, y):
                if F(s) == t:
                    sc[i] = sc[i] + rr[j] * gi
                    sG = sG - rr[j] * gi * v
                    gi = gi * gams[j]
    sW = [rr[j] * T[j] for j in range(len(T))]
    lhs = oc.msm(list(Cs) + [srs.G1] + list(proof.Ws), sc + [sG] + sW)
    rhs = oc.msm(list(proof.Ws), rr)
    return oc.pairing_check([(lhs, srs.g2[0]), (-rhs, srs.g2[1])])


# ---------------------------------------------------------------------------
# SHPLONK [BDFG20]: nhiều đa thức, mỗi đa thức một tập điểm riêng, 2 phần tử G1
# ---------------------------------------------------------------------------

ShplonkProof = namedtuple("ShplonkProof", "ys W Wp")


def _shplonk_challenge_z(tr, T):
    while True:                          # z phải nằm ngoài T (xác suất lặp ≈ |T|/r)
        z = tr.challenge(b"z")
        if z not in T:
            return z


def open_shplonk(srs, fs, Cs, Ss, tr):
    """Prover SHPLONK.

    Vòng 1: γ; f(X) = Σ γ^(i-1) Z_{T\\S_i}(X)(f_i(X) - r_i(X)); h = f/Z_T; W = [h(τ)]_1.
    Vòng 2: z; L(X) = Σ γ^(i-1) Z_{T\\S_i}(z)(f_i(X) - r_i(z)) - Z_T(z)·h(X);
            W' = [L(X)/(X - z)]_1  (L(z) = 0).
    """
    fs = [_poly(f) for f in fs]
    ys = [[f(F(s)) for s in S] for f, S in zip(fs, Ss)]
    _absorb_claims(tr, Cs, Ss, ys)
    T = _distinct(Ss)
    ZT = vanishing(T)
    rs = [lagrange_interpolate(F, S, y) for S, y in zip(Ss, ys)]
    Zc = [vanishing([t for t in T if t not in [F(s) for s in S]]) for S in Ss]
    g = _powers(tr.challenge(b"gamma"), len(fs))
    fsum = Poly.zero(F)
    for gi, Zi, f, r in zip(g, Zc, fs, rs):
        fsum = fsum + Zi * (f - r) * gi
    h = _exact_div(fsum, ZT)
    W = commit(srs, h)
    tr.absorb(b"W", W)
    z = _shplonk_challenge_z(tr, T)
    L = h * (-ZT(z))
    for gi, Zi, f, r in zip(g, Zc, fs, rs):
        L = L + (f - r(z)) * (gi * Zi(z))
    assert L(z) == 0
    Wp = commit(srs, _exact_div(L, Poly(F, [-z, 1])))
    return ShplonkProof(ys, W, Wp)


def verify_shplonk(srs, Cs, Ss, proof, tr):
    """Verifier SHPLONK: F = Σ γ^(i-1) Z_{T\\S_i}(z)(C_i - [r_i(z)]_1) - Z_T(z)·W,
    kiểm e(F + z·W', [1]_2) · e(-W', [τ]_2) = 1 — hai pairing, một MSM cỡ k + 3."""
    _absorb_claims(tr, Cs, Ss, proof.ys)
    T = _distinct(Ss)
    g = _powers(tr.challenge(b"gamma"), len(Cs))
    tr.absorb(b"W", proof.W)
    z = _shplonk_challenge_z(tr, T)
    sc, sG = [], F(0)
    for gi, S, y in zip(g, Ss, proof.ys):
        Si = [F(s) for s in S]
        zi = gi * vanishing([t for t in T if t not in Si])(z)
        ri_z = lagrange_interpolate(F, Si, y)(z)
        sc.append(zi)
        sG = sG - zi * ri_z
    lhs = oc.msm(list(Cs) + [srs.G1, proof.W, proof.Wp],
                 sc + [sG, -vanishing(T)(z), z])
    return oc.pairing_check([(lhs, srs.g2[0]), (-proof.Wp, srs.g2[1])])


# ---------------------------------------------------------------------------
# Hiding KZG (PolyCommit_Ped): C = [φ(τ)]G + [φ̂(τ)]H
# ---------------------------------------------------------------------------

def hash_to_g1(tag=b"snarklib/kzg/H", E=None):
    """Điểm H "không ai biết logarit": x = SHA256(tag‖đếm) mod q, y = căn nhỏ hơn,
    H = 12·(x, y). Trả về (H, số lần thử).

    HẠN CHẾ TRÊN TOY: nhóm chỉ có 449 phần tử nên ai cũng tìm được log_G(H)
    bằng vét cạn (xem dlog_bruteforce); tính che giấu chỉ có ý nghĩa hình thức.
    """
    if E is None:
        E = toy_generator()[0]
    ctr = 0
    while True:
        x = int.from_bytes(hashlib.sha256(tag + ctr.to_bytes(4, "big")).digest(), "big") % E.q
        ctr += 1
        v = E.rhs(x)
        if v == 0 or v.legendre() != 1:
            continue
        H = COFACTOR_TOY * E.point(x, v.sqrt())
        if not H.is_infinity():
            return H, ctr


def dlog_bruteforce(G, H, order=R):
    """Tìm β với β·G = H bằng vét cạn (chỉ khả thi vì nhóm toy có 449 phần tử)."""
    P = G.E.infinity()
    for b in range(order):
        if P == H:
            return b
        P = P + G
    return None


def hiding_setup(degree, tau, G1=None, G2=None, tag=b"snarklib/kzg/H"):
    """SRS cho hiding KZG: thêm [τ^i]H với H = hash_to_g1(tag)."""
    H, _ = hash_to_g1(tag)
    return setup(degree, tau, 1, G1, G2, H)


def commit_hiding(srs, f, fhat):
    """C = [f(τ)]_1 + [f̂(τ)]H, f̂ là đa thức ngẫu nhiên cùng bậc."""
    return oc.add(commit(srs, f), commit(srs, fhat, bases=srs.h1))


def open_hiding(srs, f, fhat, z):
    """Trả về (y, ŷ, W) với W = [ψ(τ)]_1 + [ψ̂(τ)]H."""
    f, fhat, z = _poly(f), _poly(fhat), F(z)
    y, yh = f(z), fhat(z)
    psi = _exact_div(f - y, Poly(F, [-z, 1]))
    psih = _exact_div(fhat - yh, Poly(F, [-z, 1]))
    return y, yh, oc.add(commit(srs, psi), commit(srs, psih, bases=srs.h1))


def verify_hiding(srs, C, z, y, yh, W):
    """Kiểm e(C - [y]G - [ŷ]H + z·W, [1]_2) · e(-W, [τ]_2) = 1."""
    lhs = oc.msm([C, srs.G1, srs.h1[0], W], [1, -F(y), -F(yh), F(z)])
    return oc.pairing_check([(lhs, srs.g2[0]), (-W, srs.g2[1])])


def fresh_transcript(label=b"ch03/kzg"):
    """Transcript mới; prover và verifier mỗi bên tạo một bản với cùng nhãn."""
    return Transcript(label, R)
