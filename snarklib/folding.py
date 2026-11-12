"""Chương 9 -- đệ quy, chu trình đường cong, số học không gốc và folding.

Bốn nhóm công cụ (mã in sách ở ch09/ chỉ gọi các hàm dưới đây):

1. Đường cong nhỏ, tìm chu trình và 2-chain (mục 9.2a)
    count_points(p, a, b)          #E(F_p) cho y^2 = x^3 + a x + b, đếm bằng ký hiệu Legendre
    ec_add, ec_mul, point_order    số học affine trên số nguyên mod p (None = điểm vô cực)
    embedding_degree(p, r)         k nhỏ nhất với r | p^k - 1
    find_cycles(lo, hi, b_max)     chu trình cỡ 2: E1/F_q1: y^2 = x^3 + b1, #E1 = q2 nguyên tố,
                                   E2/F_q2: y^2 = x^3 + b2, #E2 = q1
    chain_outer_supersingular(q)   đường ngoài y^2 = x^3 + x trên F_P, P = h·q - 1, có nhóm con cấp q
    chain_inner_over(p)            đường trong trên F_p có nhóm con nguyên tố lớn, bậc nhúng nhỏ

2. Số học không gốc (mục 9.2a, eq:ch9-nonnative)
    nonnative_mul_circuit(q, F, w) mạch DSL kiểm c ≡ a·b (mod q) trong trường F ≠ F_q bằng limb w bit
    nonnative_formula(l, w, wc)    công thức đếm ràng buộc theo số limb

3. Nova: folding cho relaxed R1CS (CHỈ để kiểm và lấy dữ liệu; công thức đã có ở [ZKP] Mục 13.1)
    Theo kothapalli2022nova.txt:734-742 (Định nghĩa 11) và :745-748 (u, E gộp), :798-806
    (Construction 1). Lệch bài gốc: Z = (u, x, W) thay vì (W, x, u) (thống nhất z = (1, x, w)
    của snarklib.r1cs); KHÔNG cam kết (W, E rõ), nên không succinct, không ZK.
    RelaxedR1CS(r1cs) -> .trivial(x, w), .is_satisfied(inst, wit), .cross_term(...), .fold(...)

4. HyperNova: multi-folding một LCCCS + một CCCS bằng sum-check (mục 9.2d)
    Theo HyperNova (ePrint 2023/573, CHƯA có txt -- cài theo hiểu biết): CCCS (x; w) thỏa
    Σ_i c_i ∘_{j∈S_i} M_j z = 0; LCCCS (u, x, r, v_1..v_t; w) thỏa v_j = (M_j z)~(r),
    z = (u, x, w). Gộp: sum-check cho
        g(X) = Σ_j γ^j · eq(r, X)·(M_j z_1)~(X) + γ^{t+1} · eq(β, X)·Σ_i c_i Π_{j∈S_i} (M_j z_2)~(X)
    với khẳng định Σ_j γ^j v_j; Prover gửi σ_j = (M_j z_1)~(r'), θ_j = (M_j z_2)~(r');
    Verifier kiểm s_l(r'_l) = e_1 Σ γ^j σ_j + γ^{t+1} e_2 Σ c_i Π θ_j, rút ρ và lấy
    u' = u_1 + ρ, x' = x_1 + ρ x_2, v'_j = σ_j + ρ θ_j, r' (và w' = w_1 + ρ w_2).
    Đơn giản hóa: bỏ cam kết C (witness rõ), không ZK, μ = ν = 1.
    CCSPadded(ccs) -> linearize(x, w, tr), hypernova_prove(...), hypernova_verify(...),
                      lcccs_is_satisfied(...)
"""

from .field import GF, is_prime, factorize
from .mle import eq, eq_table, evaluate
from .r1cs import matvec
from . import sumcheck as sc


# ==========================================================================
# 1. Đường cong nhỏ
# ==========================================================================

def _legendre_table(p):
    """chi[v] ∈ {-1, 0, 1} cho v = 0..p-1."""
    chi = [-1] * p
    chi[0] = 0
    for y in range(1, (p + 1) // 2):
        chi[y * y % p] = 1
    return chi


def count_points(p, a, b, chi=None):
    """#E(F_p) = 1 + Σ_x (1 + χ(x^3 + a x + b)), kể cả điểm vô cực."""
    chi = chi or _legendre_table(p)
    return 1 + p + sum(chi[(x * x * x + a * x + b) % p] for x in range(p))


def ec_add(P, Q, a, p):
    if P is None:
        return Q
    if Q is None:
        return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P == Q:
        lam = (3 * x1 * x1 + a) * pow(2 * y1, -1, p) % p
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)


def ec_mul(k, P, a, p):
    R = None
    while k:
        if k & 1:
            R = ec_add(R, P, a, p)
        P = ec_add(P, P, a, p)
        k >>= 1
    return R


def first_point(p, a, b, start=0):
    """Điểm affine đầu tiên (x nhỏ nhất ≥ start) trên E."""
    from .field import GF as _GF
    F = _GF(p)
    for x in range(start, p):
        rhs = F(x * x * x + a * x + b)
        if rhs == 0:
            return (x, 0)
        if rhs.is_qr():
            return (x, int(rhs.sqrt()))
    return None


def point_order(P, a, p, n):
    """Cấp của P, biết n·P = O (n là #E)."""
    assert ec_mul(n, P, a, p) is None
    order = n
    for ell in factorize(n):
        while order % ell == 0 and ec_mul(order // ell, P, a, p) is None:
            order //= ell
    return order


def embedding_degree(p, r, kmax=10 ** 6):
    """k nhỏ nhất với r | p^k - 1 (cấp của p trong F_r^*)."""
    x = p % r
    if x == 0:
        return None
    for k in range(1, min(kmax, r) + 1):
        if x == 1:
            return k
        x = x * p % r
    return None


def find_cycles(lo, hi, b_max=10, limit=5):
    """Chu trình cỡ 2 dạng y^2 = x^3 + b (j = 0, như Pasta/Tweedle).

    Duyệt q1 nguyên tố trong [lo, hi), q1 ≡ 1 (mod 3) (để đường j = 0 thường,
    không siêu kỳ dị), b1 = 1..b_max; nếu q2 = #E1 nguyên tố, q2 ≠ q1, tìm b2
    với #E2(F_q2) = q1. Trả về list dict (q1 < q2 để không lặp).
    """
    out = []
    for q1 in range(lo, hi):
        if q1 % 3 != 1 or not is_prime(q1):
            continue
        chi1 = _legendre_table(q1)
        for b1 in range(1, b_max + 1):
            q2 = count_points(q1, 0, b1, chi1)
            if q2 <= q1 or q2 % 3 != 1 or not is_prime(q2):
                continue
            chi2 = _legendre_table(q2)
            for b2 in range(1, b_max + 1):
                if count_points(q2, 0, b2, chi2) == q1:
                    out.append({"q1": q1, "b1": b1, "q2": q2, "b2": b2,
                                "k1": embedding_degree(q1, q2), "k2": embedding_degree(q2, q1)})
                    break
            if out and out[-1]["q1"] == q1:
                break
        if len(out) >= limit:
            break
    return out


def chain_outer_supersingular(q, hmax=200):
    """Đường ngoài E': y^2 = x^3 + x trên F_P, P = h·q - 1 nguyên tố, P ≡ 3 (mod 4).

    Khi đó E' siêu kỳ dị, #E' = P + 1 = h·q nên E' có nhóm con cấp q: trường vô
    hướng của E' (nhóm con đó) = F_q = trường cơ sở của đường trong. Bậc nhúng
    của nhóm con cấp q là 2 (vì P ≡ -1 mod q). Đếm điểm thật để kiểm.
    """
    for h in range(1, hmax + 1):
        P = h * q - 1
        if P % 4 != 3 or not is_prime(P):
            continue
        n = count_points(P, 1, 0)
        G = first_point(P, 1, 0, start=2)
        Gq = ec_mul(h, G, 1, P)
        k = 1
        while Gq is None:            # chọn điểm khác nếu h·G = O
            G = first_point(P, 1, 0, start=G[0] + 1)
            Gq = ec_mul(h, G, 1, P)
        return {"P": P, "h": h, "n": n, "G": Gq, "order": point_order(Gq, 1, P, n),
                "k": embedding_degree(P, q)}
    return None


def chain_inner_over(p, kmax=6, a_max=None):
    """Đường trong E/F_p: y^2 = x^3 + a x + b có nhóm con nguyên tố r' lớn nhất với
    bậc nhúng ≤ kmax (thân thiện pairing kiểu thu nhỏ). Trả về dict tốt nhất."""
    chi = _legendre_table(p)
    best = None
    a_max = p if a_max is None else a_max
    for a in range(a_max):
        for b in range(1, p):
            if (4 * a ** 3 + 27 * b * b) % p == 0:
                continue
            n = count_points(p, a, b, chi)
            r = max(factorize(n))
            if r <= 3 or r == p:          # r = p: đường bất thường (anomalous), bỏ
                continue
            k = embedding_degree(p, r)
            if k <= kmax and (best is None or (r, -k) > (best["r"], -best["k"])):
                best = {"a": a, "b": b, "n": n, "r": r, "k": k}
    return best


# ==========================================================================
# 2. Số học không gốc
# ==========================================================================

# Hai trường của chu trình Pasta (halo2curves-0.9.0 src/pasta/fp.rs:20, fq.rs:20;
# pallas.rs:87-92: Pallas định nghĩa trên F_p có cấp q, Vesta ngược lại)
PALLAS_P = 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001
PALLAS_Q = 0x40000000000000000000000000000000224698fc0994a8dd8c46eb2100000001


def probable_prime(n, bases=(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)):
    """Miller–Rabin với các cơ sở cố định (cho số lớn, field.is_prime chia thử quá chậm)."""
    if n < 2:
        return False
    for b in bases:
        if n % b == 0:
            return n == b
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for b in bases:
        x = pow(b, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def big_field(p):
    """Lớp phần tử F_p cho p lớn (kiểm Miller–Rabin thay vì chia thử)."""
    from .field import FieldElement, _CACHE
    if p in _CACHE:
        return _CACHE[p]
    if not probable_prime(p):
        raise ValueError("p không nguyên tố")
    cls = type("F%dbit" % p.bit_length(), (FieldElement,), {"p": p, "__slots__": ()})
    _CACHE[p] = cls
    return cls

def _limbs(v, w, l):
    return [(v >> (w * i)) & ((1 << w) - 1) for i in range(l)]


def nonnative_params(q, p, w):
    """Số limb l, cận cột và số bit carry wc cho phép nhân không gốc limb w bit.

    Cột j của a(X)·b(X) - q(X)·k(X) - c(X) có trị tuyệt đối ≤ D = 2·l·(2^w-1)^2 + 2^w;
    carry |κ_j| ≤ (D + |κ_{j-1}|)/2^w; hệ chỉ đúng (không quấn quanh F_p) nếu
    hệ số tích ≤ l·(2^w-1)^2 < p và |cột + κ_{j-1} - 2^w κ_j| < p.
    """
    l = -(-q.bit_length() // w)
    m = (1 << w) - 1
    D = 2 * l * m * m + m
    kap = 0
    for _ in range(2 * l):
        kap = (D + kap) // (1 << w) + 1
    wc = (2 * kap).bit_length()
    bound = D + kap + (1 << w) * kap
    ok = l * m * m < p and bound < p
    return {"l": l, "w": w, "wc": wc, "kappa": kap, "col_bound": bound, "ok": ok}


def nonnative_formula(l, w, wc, range_ab=True):
    """Số ràng buộc của nonnative_mul_circuit:
    limb a, b, k, c mỗi limb (w + 1) [w ràng buộc bit + 1 tổ hợp]  -> 4l(w+1) (2l(w+1) nếu a, b đã kiểm)
    tích a(X)b(X) tại 2l-1 điểm                                      -> 2l - 1
    phương trình cột                                                 -> 2l - 1
    carry 2l-2 giá trị, mỗi giá trị (wc + 1)                          -> (2l-2)(wc+1)
    """
    rng = (4 if range_ab else 2) * l * (w + 1)
    return rng + 2 * (2 * l - 1) + (2 * l - 2) * (wc + 1)


def nonnative_mul_circuit(q, F, w, range_ab=True):
    """Mạch DSL (snarklib.dsl) kiểm a·b = q·k + c (số nguyên) với limb w bit, trên trường F.

    Biến: limb a_i, b_i, k_i, c_i (bí mật), hệ số tích p_j (2l-1), carry κ_j (2l-2,
    lưu dạng κ_j + 2^(wc-1) để không âm). Chưa kiểm c < q (biểu diễn chưa rút gọn,
    như các thư viện bignum thường làm) và giả sử a, b < q (đầu vào từ phép trước).
    Trả về (circuit, prm, assign) với assign(a, b) -> dict đầu vào cho solve().
    """
    from .dsl import Circuit
    prm = nonnative_params(q, F.p, w)
    l, wc = prm["l"], prm["wc"]
    off = 1 << (wc - 1)
    c = Circuit(F)
    A = [c.private("a%d" % i) for i in range(l)]
    B = [c.private("b%d" % i) for i in range(l)]
    K = [c.private("k%d" % i) for i in range(l)]
    C = [c.private("c%d" % i) for i in range(l)]
    Pc = [c.private("p%d" % j) for j in range(2 * l - 1)]
    Kap = [c.private("kap%d" % j) for j in range(2 * l - 2)]
    for grp in ((A, B, K, C) if range_ab else (K, C)):
        for v in grp:
            c.to_bits(v, w)
    for v in Kap:
        c.to_bits(v, wc)
    # p(X) = a(X)·b(X): kiểm tại X = 0..2l-2 (mỗi điểm 1 ràng buộc)
    for t in range(2 * l - 1):
        at = sum((A[i] * pow(t, i) for i in range(l)), c.const(0))
        bt = sum((B[i] * pow(t, i) for i in range(l)), c.const(0))
        pt = sum((Pc[j] * pow(t, j) for j in range(2 * l - 1)), c.const(0))
        c.enforce(at, bt, pt, "tich")
    ql = _limbs(q, w, l)
    for j in range(2 * l - 1):
        col = Pc[j] - (C[j] if j < l else 0)
        for i in range(l):
            if 0 <= j - i < l and ql[j - i]:
                col = col - K[i] * ql[j - i]
        if j > 0:
            col = col + (Kap[j - 1] - off)
        rhs = (Kap[j] - off) * (1 << w) if j < 2 * l - 2 else c.const(0)
        c.assert_eq(col, rhs, "cot")

    def assign(a, b, cc=None, kk=None):
        kk = a * b // q if kk is None else kk
        cc = a * b - q * kk if cc is None else cc
        al, bl, kl, cl = (_limbs(v, w, l) for v in (a, b, kk, cc))
        prod = [sum(al[i] * bl[j - i] for i in range(l) if 0 <= j - i < l) for j in range(2 * l - 1)]
        env = {}
        for nm, vals in (("a", al), ("b", bl), ("k", kl), ("c", cl), ("p", prod)):
            env.update({"%s%d" % (nm, i): v for i, v in enumerate(vals)})
        kap = 0
        for j in range(2 * l - 2):
            col = prod[j] - (cl[j] if j < l else 0) - sum(
                kl[i] * ql[j - i] for i in range(l) if 0 <= j - i < l) + kap
            kap = col >> w
            env["kap%d" % j] = kap + off
        return env

    return c, prm, assign


# ==========================================================================
# 3. Nova: relaxed R1CS (z = (u, x, W))
# ==========================================================================

class RelaxedR1CS:
    """Relaxed R1CS (A Z) ∘ (B Z) = u·(C Z) + E, Z = (u, x, W) (kothapalli2022nova.txt:734-740)."""

    def __init__(self, r1cs):
        self.r = r1cs
        self.F = r1cs.F

    def Z(self, inst, W):
        return [self.F(inst["u"])] + [self.F(v) for v in inst["x"]] + [self.F(v) for v in W]

    def trivial(self, x, w):
        """R1CS -> relaxed với u = 1, E = 0 (:741-742)."""
        F = self.F
        return {"u": F(1), "x": [F(v) for v in x], "E": [F(0)] * self.r.m}, [F(v) for v in w]

    def residual(self, inst, W):
        F, r = self.F, self.r
        Z = self.Z(inst, W)
        Az, Bz, Cz = (matvec(F, M, Z) for M in (r.A, r.B, r.C))
        return [a * b - inst["u"] * c - e for a, b, c, e in zip(Az, Bz, Cz, inst["E"])]

    def is_satisfied(self, inst, W):
        return all(v == 0 for v in self.residual(inst, W))

    def cross_term(self, i1, W1, i2, W2):
        """T = AZ1∘BZ2 + AZ2∘BZ1 - u1·CZ2 - u2·CZ1 (:798-800)."""
        F, r = self.F, self.r
        Z1, Z2 = self.Z(i1, W1), self.Z(i2, W2)
        A1, B1, C1 = (matvec(F, M, Z1) for M in (r.A, r.B, r.C))
        A2, B2, C2 = (matvec(F, M, Z2) for M in (r.A, r.B, r.C))
        return [a1 * b2 + a2 * b1 - i1["u"] * c2 - i2["u"] * c1
                for a1, b1, c1, a2, b2, c2 in zip(A1, B1, C1, A2, B2, C2)]

    def fold(self, i1, W1, i2, W2, tr, label=b"nova"):
        """Một bước NIFS (Fiat–Shamir): trả về (inst, W, T, r). Không cam kết -> T gửi rõ."""
        F = self.F
        T = self.cross_term(i1, W1, i2, W2)
        for nm, inst in ((b"/U1", i1), (b"/U2", i2)):
            tr.absorb(label + nm, [int(inst["u"])] + [int(v) for v in inst["x"]] +
                      [int(v) for v in inst["E"]])
        tr.absorb(label + b"/T", [int(v) for v in T])
        r = sc._challenge(F, tr, label + b"/r")
        inst = {"u": i1["u"] + r * i2["u"],
                "x": [a + r * b for a, b in zip(i1["x"], i2["x"])],
                "E": [e1 + r * t + r * r * e2 for e1, t, e2 in zip(i1["E"], T, i2["E"])]}
        W = [a + r * b for a, b in zip(W1, W2)]
        return inst, W, T, r


# ==========================================================================
# 4. HyperNova: multi-folding LCCCS + CCCS
# ==========================================================================

def _pow2(n):
    s = 0
    while (1 << s) < n:
        s += 1
    return s


class CCSPadded:
    """CCS (snarklib.ccs) đệm m, n lên lũy thừa 2; z = (u, x, w, 0...)."""

    def __init__(self, ccs):
        self.ccs = ccs
        self.F = F = ccs.F
        self.s = _pow2(ccs.m)
        self.sp = _pow2(ccs.n)
        self.m, self.n = 1 << self.s, 1 << self.sp
        self.M = [[[F(v) for v in row] + [F(0)] * (self.n - ccs.n) for row in M] +
                  [[F(0)] * self.n for _ in range(self.m - ccs.m)] for M in ccs.M]
        self.t = len(self.M)
        self.S, self.c = ccs.S, ccs.c
        self.d = ccs.degree()

    def z(self, u, x, w):
        F = self.F
        z = [F(u)] + [F(v) for v in x] + [F(v) for v in w]
        return z + [F(0)] * (self.n - len(z))

    def Mz(self, z):
        return [matvec(self.F, M, z) for M in self.M]

    def cccs_is_satisfied(self, x, w):
        return self.ccs.is_satisfied(x, w)

    def lcccs_is_satisfied(self, inst, w):
        """v_j = (M_j z)~(r) với z = (u, x, w)."""
        Mz = self.Mz(self.z(inst["u"], inst["x"], w))
        return all(evaluate(self.F, mz, inst["r"]) == v for mz, v in zip(Mz, inst["v"]))

    def linearize(self, x, w, tr, label=b"lin"):
        """CCCS thỏa -> LCCCS với u = 1, r rút từ transcript (Prover trung thực).

        Trong HyperNova bước này cũng là một sum-check (multi-folding với μ = 0);
        ở đây chỉ dùng để dựng instance tích lũy ban đầu.
        """
        F = self.F
        tr.absorb(label + b"/x", [int(v) for v in x])
        r = [sc._challenge(F, tr, label + b"/r%d" % i) for i in range(self.s)]
        Mz = self.Mz(self.z(1, x, w))
        return {"u": F(1), "x": [F(v) for v in x], "r": r,
                "v": [evaluate(F, mz, r) for mz in Mz]}

    def _g(self, gamma):
        t = self.t
        terms = [(gamma ** (j + 1), [1, 2 + j]) for j in range(t)]
        gt1 = gamma ** (t + 1)
        terms += [(gt1 * ci, [0] + [2 + t + j for j in Si]) for ci, Si in zip(self.c, self.S)]
        return sc.sum_of_products(terms)

    def _challenges(self, tr, label, l1, x2):
        F = self.F
        tr.absorb(label + b"/L", [int(l1["u"])] + [int(v) for v in l1["x"] + l1["r"] + l1["v"]])
        tr.absorb(label + b"/C", [int(v) for v in x2])
        gamma = sc._challenge(F, tr, label + b"/gamma")
        beta = [sc._challenge(F, tr, label + b"/beta%d" % i) for i in range(self.s)]
        return gamma, beta

    def hypernova_prove(self, l1, w1, x2, w2, tr, label=b"hn"):
        """Gộp LCCCS (l1; w1) và CCCS (x2; w2). Trả về (proof, l_new, w_new)."""
        F = self.F
        gamma, beta = self._challenges(tr, label, l1, x2)
        Mz1 = self.Mz(self.z(l1["u"], l1["x"], w1))
        Mz2 = self.Mz(self.z(1, x2, w2))
        tables = [eq_table(F, beta), eq_table(F, l1["r"])] + Mz1 + Mz2
        g, deg = self._g(gamma)
        claim = sum((gamma ** (j + 1) * v for j, v in enumerate(l1["v"])), F(0))
        proof, rp, finals = sc.prove(F, tables, g, deg, tr, claim=claim, label=label + b"/sc")
        sigma, theta = finals[2:2 + self.t], finals[2 + self.t:]
        tr.absorb(label + b"/st", [int(v) for v in sigma + theta])
        rho = sc._challenge(F, tr, label + b"/rho")
        l_new = {"u": l1["u"] + rho, "x": [a + rho * b for a, b in zip(l1["x"], x2)],
                 "r": rp, "v": [s + rho * th for s, th in zip(sigma, theta)]}
        w_new = [F(a) + rho * b for a, b in zip(w1, w2)]
        return {"sc": proof, "sigma": sigma, "theta": theta}, l_new, w_new

    def hypernova_verify(self, l1, x2, proof, tr, label=b"hn"):
        """Verifier gộp: trả về (ok, l_new). Chỉ dùng instance, không dùng witness."""
        F = self.F
        gamma, beta = self._challenges(tr, label, l1, x2)
        _, deg = self._g(gamma)
        claim = sum((gamma ** (j + 1) * v for j, v in enumerate(l1["v"])), F(0))
        ok, rp, final = sc.verify(F, self.s, deg, proof["sc"], tr, label=label + b"/sc", claim=claim)
        if not ok:
            return False, None
        sigma, theta = proof["sigma"], proof["theta"]
        e1, e2 = eq(F, l1["r"], rp), eq(F, beta, rp)
        rhs = e1 * sum((gamma ** (j + 1) * s for j, s in enumerate(sigma)), F(0))
        q = F(0)
        for ci, Si in zip(self.c, self.S):
            prod = ci
            for j in Si:
                prod = prod * theta[j]
            q = q + prod
        rhs = rhs + gamma ** (self.t + 1) * e2 * q
        if rhs != final:
            return False, None
        tr.absorb(label + b"/st", [int(v) for v in list(sigma) + list(theta)])
        rho = sc._challenge(F, tr, label + b"/rho")
        return True, {"u": l1["u"] + rho, "x": [a + rho * b for a, b in zip(l1["x"], x2)],
                      "r": rp, "v": [s + rho * th for s, th in zip(sigma, theta)]}
