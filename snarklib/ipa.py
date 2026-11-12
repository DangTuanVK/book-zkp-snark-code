"""Inner Product Argument (IPA) kiểu Bulletproofs và PCS dựng từ IPA.

Ký hiệu cộng trên nhóm con bậc r = 449 của đường cong toy (chỉ dùng phép cộng
điểm, không dùng pairing). Véc-tơ cơ sở lấy từ snarklib.pedersen.

1) IPA tổng quát (Bünz et al. 2018, Protocol 1 + Protocol 2): chứng minh biết
   a, b với P = <a, G> + <b, H> và c = <a, b>.
   Protocol 1:  x ← thách thức,  P' = P + (x·c)·U,  U' = x·U.
   Protocol 2 (mỗi vòng, n' = n/2, lo = nửa đầu, hi = nửa sau):
       c_L = <a_lo, b_hi>,   c_R = <a_hi, b_lo>
       L = <a_lo, G_hi> + <b_hi, H_lo> + c_L·U
       R = <a_hi, G_lo> + <b_lo, H_hi> + c_R·U
       x ← thách thức (Fiat–Shamir, sau khi hấp thụ L, R)
       G' = x^-1·G_lo + x·G_hi,     H' = x·H_lo + x^-1·H_hi
       P' = x^2·L + P + x^-2·R
       a' = x·a_lo + x^-1·a_hi,     b' = x^-1·b_lo + x·b_hi
   Khi n = 1: gửi a, b; kiểm P = a·G + b·H + (a·b)·U.
   (Đối chiếu bunz2018bulletproofs.txt dòng 813–862 và 755–766.)

2) IPA với b CÔNG KHAI (biến thể của Halo, bowe2020recursive.txt dòng 333–398):
   không cần H; L = <a_lo, G_hi> + <a_lo, b_hi>·U, R = <a_hi, G_lo> + <a_hi, b_lo>·U.
   Verifier tự gập b. Proof = (L_j, R_j)_{j=1..k} và một vô hướng a cuối.

3) PCS từ IPA: Commit(p) = <a, G> với a = hệ số của p (không ẩn, không blinding);
   mở tại z: v = p(z) = <a, b>, b = (1, z, z^2, ..., z^{n-1}).

Verifier "một MSM" (Bünz et al. Mục 3.1): G cuối = <s, G> với
   s_i = Π_j u_j^{±1}  (+1 nếu bit thứ j của i, tính từ bit CAO, bằng 1).
"""

from .field import GF
from .curve import R_TOY
from .pedersen import E, PedersenParams, msm, smul, padd

F = GF(R_TOY)


# ----------------------------------------------------------------- tiện ích
def inner(a, b):
    return sum((F(x) * F(y) for x, y in zip(a, b)), F(0))


def powers(z, n):
    z = F(z)
    out, t = [], F(1)
    for _ in range(n):
        out.append(t)
        t = t * z
    return out


def _log2(n):
    if n == 0 or n & (n - 1):
        raise ValueError("n phải là lũy thừa của 2")
    return n.bit_length() - 1


def challenge_nonzero(t, label):
    """Thách thức khác 0 (cần nghịch đảo); rút lại với cùng nhãn nếu ra 0."""
    while True:
        x = t.challenge(label)
        if x != 0:
            return x


def s_vector(us, cnt=None):
    """s_i = Π_j u_j^{+1 nếu bit j (từ cao) của i là 1, ngược lại -1}; n nhân trường."""
    s = [F(1)]
    for u in us:                        # mỗi vòng sau ứng với một bit thấp hơn
        ui = u.inverse()
        s = [x * w for x in s for w in (ui, u)]
        if cnt is not None:
            cnt.fmul += len(s)
    return s


def b_final_powers(us, z):
    """<s, (1, z, ..., z^{n-1})> = Π_j (u_j^{-1} + u_j·z^{2^{k-j}}) — O(log n)."""
    k = len(us)
    acc = F(1)
    for j, u in enumerate(us):
        acc = acc * (u.inverse() + u * F(z) ** (1 << (k - 1 - j)))
    return acc


class VCount:
    """Đếm phép toán Verifier: nhân vô hướng, cộng điểm, nhân trường (ước lượng)."""

    def __init__(self):
        self.smul = self.add = self.fmul = 0

    def as_dict(self):
        return {"smul": self.smul, "add": self.add, "fmul": self.fmul}


# ------------------------------------------------- một vòng (b công khai)
def round_messages(G, a, b, U):
    """Thông điệp một vòng: (L, R, c_L, c_R) với b công khai."""
    h = len(a) // 2
    cL = inner(a[:h], b[h:])
    cR = inner(a[h:], b[:h])
    L = padd(msm(a[:h], G[h:]), smul(U, cL))
    R = padd(msm(a[h:], G[:h]), smul(U, cR))
    return L, R, cL, cR


def fold_vectors(G, a, b, x):
    """Gập đôi: G' = x^-1 G_lo + x G_hi; a' = x a_lo + x^-1 a_hi; b' = x^-1 b_lo + x b_hi."""
    h, xi = len(a) // 2, x.inverse()
    G2 = [smul(G[i], xi) + smul(G[h + i], x) for i in range(h)]
    a2 = [x * a[i] + xi * a[h + i] for i in range(h)]
    b2 = [xi * b[i] + x * b[h + i] for i in range(h)]
    return G2, a2, b2


def fold_commitment(P, L, R, x, cnt=None):
    """P' = x^2·L + P + x^-2·R."""
    x2 = x * x
    return padd(padd(smul(L, x2, cnt), P, cnt), smul(R, x2.inverse(), cnt), cnt)


# ------------------------------------------------- IPA với b công khai
def prove_public_b(G, U, a, b, transcript):
    """Chứng minh P = <a, G> + <a, b>·U với b công khai. Trả về proof (dict)."""
    a = [F(v) for v in a]
    b = [F(v) for v in b]
    G = list(G)
    _log2(len(a))
    Ls, Rs = [], []
    while len(a) > 1:
        L, R, _, _ = round_messages(G, a, b, U)
        transcript.absorb(b"L", L)
        transcript.absorb(b"R", R)
        x = challenge_nonzero(transcript, b"u")
        G, a, b = fold_vectors(G, a, b, x)
        Ls.append(L)
        Rs.append(R)
    return {"L": Ls, "R": Rs, "a": a[0]}


def verify_public_b(G, U, P, b, proof, transcript, cnt=None, b_final=None):
    """Kiểm bằng một MSM: Q = P + Σ u_j^2 L_j + Σ u_j^-2 R_j  =?  a·(<s,G>) + (a·b_f)·U.

    b_final: hàm us -> <s, b> nếu b có cấu trúc (vd. lũy thừa của z); mặc định O(n).
    """
    n = len(G)
    k = _log2(n)
    if len(proof["L"]) != k or len(proof["R"]) != k:
        return False
    us = []
    Q = P
    for L, R in zip(proof["L"], proof["R"]):
        transcript.absorb(b"L", L)
        transcript.absorb(b"R", R)
        u = challenge_nonzero(transcript, b"u")
        us.append(u)
        Q = fold_commitment(Q, L, R, u, cnt)
    s = s_vector(us, cnt)
    G_f = msm(s, G, cnt)
    if b_final is not None:
        b_f = b_final(us)
    else:
        b_f = inner(s, b)
        if cnt is not None:
            cnt.fmul += n
    a = F(proof["a"])
    rhs = padd(smul(G_f, a, cnt), smul(U, a * b_f, cnt), cnt)
    return Q == rhs


# ------------------------------------------------- PCS từ IPA
def setup(n, label=b"snarklib/ipa"):
    return PedersenParams(n, label)


def commit(pp, coeffs):
    """Cam kết đa thức bậc < n: P = <a, G> (không blinding)."""
    return msm([F(c) for c in coeffs], pp.G[:len(coeffs)])


def _bind(transcript, pp, P, z, v):
    transcript.absorb(b"P", P)
    transcript.absorb(b"z", F(z))
    transcript.absorb(b"v", F(v))
    xu = challenge_nonzero(transcript, b"x_u")       # Protocol 1: U' = x_u·U
    return smul(pp.U, xu)


def open(pp, coeffs, z, transcript):
    """Mở tại z. Trả về (v, proof)."""
    a = [F(c) for c in coeffs]
    n = len(a)
    b = powers(z, n)
    v = inner(a, b)
    P = commit(pp, a)
    U1 = _bind(transcript, pp, P, z, v)
    return v, prove_public_b(pp.G[:n], U1, a, b, transcript)


def verify(pp, P, z, v, proof, transcript, n=None, cnt=None):
    """Kiểm proof mở P tại z ra v. n = số hệ số (mặc định pp.n)."""
    n = n or pp.n
    if cnt is not None:
        cnt.smul += 1          # x_u·U
        cnt.smul += 1          # v·U'
        cnt.add += 1
    U1 = _bind(transcript, pp, P, z, v)
    P1 = padd(P, smul(U1, v))                       # P' = P + v·U'
    return verify_public_b(pp.G[:n], U1, P1, powers(z, n), proof, transcript, cnt,
                           b_final=lambda us: b_final_powers(us, z))


# ------------------------------------------------- IPA tổng quát (Bulletproofs)
def bp_prove(G, H, U, a, b, transcript):
    """Bulletproofs Protocol 1 + 2 cho P = <a,G> + <b,H>, c = <a,b>. Trả về proof."""
    a = [F(v) for v in a]
    b = [F(v) for v in b]
    G, H = list(G), list(H)
    _log2(len(a))
    c = inner(a, b)
    P = padd(msm(a, G), msm(b, H))
    transcript.absorb(b"P", P)
    transcript.absorb(b"c", c)
    x = challenge_nonzero(transcript, b"x")
    U = smul(U, x)                                     # U' = x·U
    Ls, Rs = [], []
    while len(a) > 1:
        h = len(a) // 2
        cL, cR = inner(a[:h], b[h:]), inner(a[h:], b[:h])
        L = padd(padd(msm(a[:h], G[h:]), msm(b[h:], H[:h])), smul(U, cL))
        R = padd(padd(msm(a[h:], G[:h]), msm(b[:h], H[h:])), smul(U, cR))
        transcript.absorb(b"L", L)
        transcript.absorb(b"R", R)
        x = challenge_nonzero(transcript, b"x")
        xi = x.inverse()
        G = [smul(G[i], xi) + smul(G[h + i], x) for i in range(h)]
        H = [smul(H[i], x) + smul(H[h + i], xi) for i in range(h)]
        a = [x * a[i] + xi * a[h + i] for i in range(h)]
        b = [xi * b[i] + x * b[h + i] for i in range(h)]
        Ls.append(L)
        Rs.append(R)
    return {"L": Ls, "R": Rs, "a": a[0], "b": b[0]}


def bp_verify(G, H, U, P, c, proof, transcript, method="msm", cnt=None):
    """Kiểm Bulletproofs IPA. method = 'fold' (gập từng vòng, 4n nhân vô hướng)
    hoặc 'msm' (một MSM cỡ 2n + 2 log n + 1, Bünz et al. Mục 3.1)."""
    n = len(G)
    k = _log2(n)
    if len(proof["L"]) != k or len(proof["R"]) != k:
        return False
    transcript.absorb(b"P", P)
    transcript.absorb(b"c", F(c))
    x = challenge_nonzero(transcript, b"x")
    U = smul(U, x, cnt)
    P = padd(P, smul(U, c, cnt), cnt)                 # P' = P + (x·c)·U
    G, H = list(G), list(H)
    us = []
    for L, R in zip(proof["L"], proof["R"]):
        transcript.absorb(b"L", L)
        transcript.absorb(b"R", R)
        u = challenge_nonzero(transcript, b"x")
        us.append(u)
        P = fold_commitment(P, L, R, u, cnt)
        if method == "fold":
            ui, h = u.inverse(), len(G) // 2
            G = [padd(smul(G[i], ui, cnt), smul(G[h + i], u, cnt), cnt) for i in range(h)]
            H = [padd(smul(H[i], u, cnt), smul(H[h + i], ui, cnt), cnt) for i in range(h)]
    a, b = F(proof["a"]), F(proof["b"])
    if method == "fold":
        Gf, Hf = G[0], H[0]
    else:
        s = s_vector(us, cnt)
        Gf = msm(s, G, cnt)
        Hf = msm([v.inverse() for v in s], H, cnt)     # h cuối dùng 1/s_i
        if cnt is not None:
            cnt.fmul += n
    rhs = padd(padd(smul(Gf, a, cnt), smul(Hf, b, cnt), cnt), smul(U, a * b, cnt), cnt)
    return P == rhs


# ------------------------------------------------- kích thước
def proof_size_pcs(n):
    """Proof PCS-IPA: 2·log2 n điểm + 1 vô hướng (a cuối)."""
    return {"group": 2 * _log2(n), "field": 1}


def proof_size_bp(n):
    """Proof Bulletproofs IPA (Protocol 2): 2·log2 n điểm + 2 vô hướng."""
    return {"group": 2 * _log2(n), "field": 2}
