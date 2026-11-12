"""Giao thức sum-check (Lund–Fortnow–Karloff–Nisan) và hai PIOP dựng trên nó.

Bài toán: Prover thuyết phục Verifier rằng
    H = Σ_{b ∈ {0,1}^l}  g( T_1~(b), T_2~(b), ..., T_K~(b) )
trong đó T_k~ là MLE của bảng T_k (2^l giá trị, thứ tự big-endian như
snarklib.mle: x_1 là bit cao, được cố định TRƯỚC), còn g là một đa thức K biến
tùy ý có bậc tổng ≤ d theo các đối số (vd. g(u, v) = u·v, g(e, a, b, c) =
e·(a·b - c)). Khi đó biểu thức dưới dấu tổng có bậc ≤ d theo từng biến x_j.

Giao thức (đối chiếu thaler2022proofs.txt:1578-1610, "Description of
Sum-Check Protocol"): ở vòng j, Prover gửi đa thức một biến
    s_j(X) = Σ_{x_{j+1..l} ∈ {0,1}} g(..)(r_1, ..., r_{j-1}, X, x_{j+1}, ..., x_l)
dưới dạng d + 1 giá trị s_j(0), s_j(1), ..., s_j(d); Verifier kiểm
s_j(0) + s_j(1) = s_{j-1}(r_{j-1}) (vòng 1: = H), rút r_j từ transcript. Cuối
cùng Verifier cần MỘT lần truy vấn oracle: g(T_1~(r), ..., T_K~(r)) = s_l(r_l).
Soundness: lỗi ≤ l·d/|F| (thaler2022proofs.txt:1713-1720, Proposition 4.1).

Prover thời gian tuyến tính bằng bảng gập đôi (thaler2022proofs.txt:2347-2357,
Lemma 4.3; :2361 Lemma 4.4): ở vòng j, với mỗi bảng (đã cố định x_1..x_{j-1})
    T(t, x') = T(0, x') + t·(T(1, x') - T(0, x'))     (t = 0..d, tính bằng phép CỘNG)
và sau khi có r_j:  T'[x'] = T[x'] + r_j·(T[x' + half] - T[x'])   (1 phép nhân/ô).
Tổng chi phí ≈ Σ_j 2^{l-j}·(K + (d+1)·cost(g)) = O(2^l) phép nhân.

Proof gửi đủ d + 1 giá trị mỗi vòng (l·(d+1) phần tử trường) cho dễ đọc; bản
tối ưu bỏ s_j(1) vì Verifier tự suy ra từ s_j(1) = khẳng định - s_j(0).

API chính
    prove(F, tables, g, degree, transcript, claim=None, label=b"sc")
        -> (proof, r, finals)   proof = {"claim", "polys": [[s_j(0..d)], ...]}
                                finals[k] = T_k~(r)  (Prover dùng để trả lời oracle)
    verify(F, l, degree, proof, transcript, label=b"sc", claim=None)
        -> (ok, r, final)       ok: các kiểm tra vòng; final = s_l(r_l), người gọi
                                PHẢI kiểm thêm g(oracle(r)) == final.
    sum_of_products(terms)      -> (g, degree): g = Σ c_i · Π_{k ∈ S_i} u_k
    zerocheck_prove/verify      f(x) = 0 trên mọi x ∈ {0,1}^l  (HyperPlonk §3.2)
    productcheck_prove/verify   Π_x f1(x)/f2(x) = s           (HyperPlonk §3.3)
    plonkish_tables, hyperplonk_prove/verify
                                PIOP HyperPlonk tối giản cho mạch Plonkish (cổng +
                                hoán vị), oracle là chính các bảng (chưa có PCS)
    counting_field(F)           lớp con của F đếm phép nhân (để lập bảng chi phí)
"""

import copy

from .mle import eq_table, eq, num_vars


# --------------------------------------------------------------------------
# Trường đếm phép nhân
# --------------------------------------------------------------------------

def counting_field(F):
    """Lớp con của F đếm mọi phép nhân (*) và nghịch đảo trong cls.counter.

    Phép cộng/trừ và lũy thừa (pow) không đếm. Dùng: CF = counting_field(F);
    đưa mọi bảng về CF(...), chạy thuật toán, đọc CF.counter["mul"].
    """
    box = {"mul": 0, "inv": 0}

    def __mul__(self, o):
        box["mul"] += 1
        return F.__mul__(self, o)

    def inverse(self):
        box["inv"] += 1
        return F.inverse(self)

    cls = type("Dem" + F.__name__, (F,),
               {"__slots__": (), "__mul__": __mul__, "__rmul__": __mul__,
                "inverse": inverse, "counter": box})
    return cls


def mul_count(F):
    """Số phép nhân đã đếm (0 nếu F không phải trường đếm)."""
    return getattr(F, "counter", {}).get("mul", 0)


# --------------------------------------------------------------------------
# Tiện ích
# --------------------------------------------------------------------------

def interpolate_eval(F, evals, x):
    """Giá trị tại x của đa thức bậc ≤ d cho bởi evals = [p(0), ..., p(d)] (Lagrange)."""
    d = len(evals) - 1
    x = F(x)
    acc = F(0)
    for i in range(d + 1):
        num, den = F(1), F(1)
        for j in range(d + 1):
            if j != i:
                num = num * (x - j)
                den = den * (i - j)
        acc = acc + F(evals[i]) * num / den
    return acc


def sum_of_products(terms):
    """terms = [(c_i, [k, k', ...]), ...] -> (g, degree) với g(u) = Σ c_i Π u_k.

    Dạng này đủ cho R1CS (e·(a·b - c)), CCS/HyperNova (Σ c_i Π_{j∈S_i} (M_j z)~),
    cổng PLONK... Bậc = max |S_i|.
    """
    terms = [(c, list(S)) for c, S in terms]

    def g(u):
        acc = None
        for c, S in terms:
            p = None
            for k in S:
                p = u[k] if p is None else p * u[k]
            if p is None:
                p = c
            elif c != 1:
                p = p * c
            acc = p if acc is None else acc + p
        return acc

    return g, max(len(S) for _, S in terms)


def soundness_bound(l, d, p):
    """Cận lỗi l·d/p của sum-check (Proposition 4.1 trong [Thaler])."""
    return l * d / p


def _challenge(F, transcript, label):
    return F(int(transcript.challenge(label)))


def _absorb_header(F, transcript, label, l, degree, claim):
    transcript.absorb(label + b"/hdr", [l, degree])
    transcript.absorb(label + b"/claim", int(claim))


# --------------------------------------------------------------------------
# Sum-check
# --------------------------------------------------------------------------

def prove(F, tables, g, degree, transcript, claim=None, label=b"sc"):
    """Prover sum-check cho Σ_b g(T_1~(b), ..., T_K~(b)).

    tables: K bảng cùng độ dài 2^l (sẽ được sao chép, không sửa bảng gốc).
    g: hàm nhận list K phần tử trường; degree: bậc theo mỗi biến (≥ 1).
    claim: nếu None, Prover tự tính tổng thật.
    """
    T = [[F(v) for v in t] for t in tables]
    l = num_vars(T[0])
    if any(len(t) != len(T[0]) for t in T):
        raise ValueError("các bảng phải cùng độ dài")
    if claim is None:
        claim = F(0)
        for i in range(len(T[0])):
            claim = claim + g([t[i] for t in T])
    claim = F(claim)
    _absorb_header(F, transcript, label, l, degree, claim)
    polys, r = [], []
    for j in range(l):
        half = len(T[0]) // 2
        s = [F(0)] * (degree + 1)
        for i in range(half):
            cur = [t[i] for t in T]                       # giá trị tại X = 0
            diff = [t[i + half] - t[i] for t in T]        # độ dốc theo X
            s[0] = s[0] + g(cur)
            for x in range(1, degree + 1):
                cur = [c + dd for c, dd in zip(cur, diff)]   # chỉ phép cộng
                s[x] = s[x] + g(cur)
        polys.append(s)
        transcript.absorb(label + b"/s%d" % j, [int(v) for v in s])
        rj = _challenge(F, transcript, label + b"/r%d" % j)
        r.append(rj)
        T = [[t[i] + rj * (t[i + half] - t[i]) for i in range(half)] for t in T]   # gập đôi
    finals = [t[0] for t in T]
    return {"claim": claim, "polys": polys}, r, finals


def verify(F, l, degree, proof, transcript, label=b"sc", claim=None):
    """Verifier sum-check. Trả về (ok, r, final); người gọi kiểm g(oracle(r)) == final.

    claim (nếu cho) phải trùng proof["claim"]; ví dụ zerocheck đòi claim = 0.
    """
    H = F(proof["claim"])
    if claim is not None and H != F(claim):
        return False, None, None
    polys = proof["polys"]
    if len(polys) != l or any(len(s) != degree + 1 for s in polys):
        return False, None, None
    _absorb_header(F, transcript, label, l, degree, H)
    cur, r = H, []
    for j, s in enumerate(polys):
        s = [F(v) for v in s]
        if s[0] + s[1] != cur:
            return False, None, None
        transcript.absorb(label + b"/s%d" % j, [int(v) for v in s])
        rj = _challenge(F, transcript, label + b"/r%d" % j)
        r.append(rj)
        cur = interpolate_eval(F, s, rj)
    return True, r, cur


def verify_with_oracle(F, l, degree, proof, transcript, g, oracle, label=b"sc", claim=None):
    """verify + kiểm cuối g(oracle(r)) == s_l(r_l); oracle(r) trả list T_k~(r)."""
    ok, r, final = verify(F, l, degree, proof, transcript, label, claim)
    return bool(ok and g(oracle(r)) == final), r


def prover_cost(l, K, degree, g_muls):
    """Số phép nhân lý thuyết của prove(): Σ_j 2^{l-j}·((d+1)·g_muls + K)."""
    return sum((1 << (l - j - 1)) * ((degree + 1) * g_muls + K) for j in range(l))


# --------------------------------------------------------------------------
# ZeroCheck (HyperPlonk, chen2023hyperplonk.txt:948-981)
# --------------------------------------------------------------------------

def zerocheck_prove(F, tables, g, degree, transcript, label=b"zc"):
    """Chứng minh g(T~(x)) = 0 với mọi x ∈ {0,1}^l bằng sum-check của eq(x, τ)·g.

    Trả về (proof, r, finals) như prove(); finals không gồm bảng eq.
    """
    l = num_vars(tables[0])
    transcript.absorb(label + b"/l", l)
    tau = [_challenge(F, transcript, label + b"/tau%d" % i) for i in range(l)]
    E = eq_table(F, tau)
    proof, r, finals = prove(F, [E] + list(tables), lambda u: u[0] * g(u[1:]),
                             degree + 1, transcript, claim=F(0), label=label)
    return proof, r, finals[1:]


def zerocheck_verify(F, l, degree, proof, transcript, label=b"zc"):
    """Trả về (ok, r, final, eq_val): người gọi kiểm eq_val·g(oracle(r)) == final."""
    transcript.absorb(label + b"/l", l)
    tau = [_challenge(F, transcript, label + b"/tau%d" % i) for i in range(l)]
    ok, r, final = verify(F, l, degree + 1, proof, transcript, label, claim=0)
    if not ok:
        return False, None, None, None
    return True, r, final, eq(F, tau, r)


# --------------------------------------------------------------------------
# ProductCheck (HyperPlonk §3.3, kỹ thuật Quarks; chen2023hyperplonk.txt:1030-1060)
# --------------------------------------------------------------------------

def product_tree(F, frac):
    """Bảng ṽ (2^{l+1} phần tử): ṽ(0,x) = frac(x), ṽ(1,x) = ṽ(x,0)·ṽ(x,1), ṽ(1..1) = 0.

    Big-endian: ṽ(0,x) ở chỉ số x, ṽ(1,x) ở 2^l + x, ṽ(x,0) ở 2x, ṽ(x,1) ở 2x + 1;
    ṽ(1,...,1,0) ở chỉ số 2^{l+1} - 2 bằng tích mọi frac(x).
    """
    n = len(frac)
    v = [F(x) for x in frac] + [F(0)] * n
    for x in range(n - 1):
        v[n + x] = v[2 * x] * v[2 * x + 1]
    return v


def _split_v(v):
    n = len(v) // 2
    V1 = v[n:]                                  # ṽ(1, x)
    VE = [v[2 * x] for x in range(n)]           # ṽ(x, 0)
    VO = [v[2 * x + 1] for x in range(n)]       # ṽ(x, 1)
    V0 = v[:n]                                  # ṽ(0, x)
    return V1, VE, VO, V0


def productcheck_prove(F, tables, f1, f2, degree, transcript, label=b"pc"):
    """Chứng minh Π_{x∈{0,1}^l} f1(T~(x))/f2(T~(x)) = s.

    f1, f2: hàm trên list giá trị các bảng, bậc ≤ degree. Prover tính ṽ, sau đó
    chạy MỘT zerocheck µ biến, gộp hai đẳng thức bằng α ngẫu nhiên
    (kỹ thuật batching, chen2023hyperplonk.txt:979-981):
        [ṽ(1,x) - ṽ(x,0)ṽ(x,1)] + α·[f2(x)·ṽ(0,x) - f1(x)] = 0.
    Bài gốc ghép hai đa thức bằng merge (µ + 1 biến, :784-786, :1054); hai cách
    tương đương về ý tưởng, cách gộp α đỡ một biến.
    Trả về (proof, r, finals, v): proof["s"] là tích = ṽ(1,..,1,0),
    finals = [ṽ(1,r), ṽ(r,0), ṽ(r,1), ṽ(0,r)] + [T_k~(r)], v là bảng ṽ.
    Lưu ý: ở đây trả lời oracle của ṽ lấy từ bảng; hệ thật cam kết ṽ bằng PCS.
    """
    T = [[F(x) for x in t] for t in tables]
    n = len(T[0])
    frac = [f1([t[i] for t in T]) / f2([t[i] for t in T]) for i in range(n)]
    v = product_tree(F, frac)
    s = v[2 * n - 2]
    transcript.absorb(label + b"/s", int(s))
    alpha = _challenge(F, transcript, label + b"/alpha")
    V1, VE, VO, V0 = _split_v(v)

    def h(u):
        return (u[0] - u[1] * u[2]) + alpha * (f2(u[4:]) * u[3] - f1(u[4:]))

    proof, r, finals = zerocheck_prove(F, [V1, VE, VO, V0] + T, h,
                                       max(2, degree + 1), transcript, label)
    proof["s"] = s
    return proof, r, finals, v


def productcheck_verify(F, l, degree, proof, transcript, f1, f2, oracle_v, oracle_t,
                        label=b"pc"):
    """oracle_v(point l+1 chiều) -> ṽ(point); oracle_t(r) -> [T_k~(r)].

    Trả về (ok, s). Kiểm: zerocheck đúng, giá trị cuối khớp, ṽ(1,..,1,0) = s.
    """
    s = F(proof["s"])
    transcript.absorb(label + b"/s", int(s))
    alpha = _challenge(F, transcript, label + b"/alpha")
    ok, r, final, eqv = zerocheck_verify(F, l, max(2, degree + 1), proof, transcript, label)
    if not ok:
        return False, s
    one, zero = F(1), F(0)
    v1 = oracle_v([one] + r)
    ve = oracle_v(r + [zero])
    vo = oracle_v(r + [one])
    v0 = oracle_v([zero] + r)
    t = oracle_t(r)
    hval = (v1 - ve * vo) + alpha * (f2(t) * v0 - f1(t))
    if eqv * hval != final:
        return False, s
    if oracle_v([one] * l + [zero]) != s:
        return False, s
    return True, s


# --------------------------------------------------------------------------
# HyperPlonk tối giản: cổng (zerocheck) + hoán vị (productcheck)
# --------------------------------------------------------------------------

GATE_KEYS = ("qL", "qR", "qO", "qM", "qCPI", "a", "b", "c")
PERM_KEYS = ("a", "b", "c", "id_a", "id_b", "id_c", "s_a", "s_b", "s_c")


def plonkish_tables(F, circ, a, b, c, x=()):
    """Bảng (độ dài 2^µ, đệm 0) cho mạch snarklib.plonkish và ba cột witness.

    qCPI = q_C + PI (gộp hằng và đầu vào công khai). id_col(x) = j·2^µ + x với j là
    chỉ số cột; s_col = σ(id) dựng từ ràng buộc sao chép (mỗi lớp tương đương
    thành một chu trình).
    """
    n = 1
    while n < circ.n:
        n *= 2
    pad = lambda col: [F(v) for v in col] + [F(0)] * (n - len(col))
    PI = circ.public_input_column(list(x))
    T = {k: pad(circ.sel[k]) for k in ("qL", "qR", "qO", "qM")}
    T["qCPI"] = pad([qc + pi for qc, pi in zip(circ.sel["qC"], PI)])
    T["a"], T["b"], T["c"] = pad(a), pad(b), pad(c)
    cid = lambda col, row: "abc".index(col) * n + row
    parent = list(range(3 * n))

    def find(u):
        while parent[u] != u:
            parent[u] = parent[parent[u]]
            u = parent[u]
        return u

    for (c1, r1), (c2, r2) in circ.copies:
        parent[find(cid(c1, r1))] = find(cid(c2, r2))
    classes = {}
    for u in range(3 * n):
        classes.setdefault(find(u), []).append(u)
    sigma = list(range(3 * n))
    for cls in classes.values():
        for i, u in enumerate(cls):
            sigma[u] = cls[(i + 1) % len(cls)]
    for j, col in enumerate("abc"):
        T["id_" + col] = [F(j * n + i) for i in range(n)]
        T["s_" + col] = [F(sigma[j * n + i]) for i in range(n)]
    return T


def _gate(u):
    qL, qR, qO, qM, qC, a, b, c = u
    return qL * a + qR * b + qO * c + qM * a * b + qC


def _perm_fns(beta, gamma):
    def f1(u):
        return (u[0] + beta * u[3] + gamma) * (u[1] + beta * u[4] + gamma) * (u[2] + beta * u[5] + gamma)

    def f2(u):
        return (u[0] + beta * u[6] + gamma) * (u[1] + beta * u[7] + gamma) * (u[2] + beta * u[8] + gamma)
    return f1, f2


def hyperplonk_prove(F, T, transcript, label=b"hp"):
    """PIOP HyperPlonk: zerocheck cổng (bậc 3) + productcheck hoán vị (s = 1)."""
    pg, rg, _ = zerocheck_prove(F, [T[k] for k in GATE_KEYS], _gate, 3, transcript, label + b"/gate")
    # Productcheck đòi f2(x) ≠ 0 trên siêu lập phương. Trên F_449 với 3·2^µ thừa số,
    # xác suất có thừa số 0 không nhỏ (≈ 3·2^µ/449), nên Prover thử lại với nonce
    # (ghi vào proof). Ở trường 256 bit xác suất này không đáng kể.
    P = [[F(v) for v in T[k]] for k in PERM_KEYS]
    for nonce in range(256):
        t2 = copy.deepcopy(transcript)
        t2.absorb(label + b"/nonce", nonce)
        beta = _challenge(F, t2, label + b"/beta")
        gamma = _challenge(F, t2, label + b"/gamma")
        f1, f2 = _perm_fns(beta, gamma)
        if all(f2([t[i] for t in P]) != 0 for i in range(len(P[0]))):
            break
    transcript.state = t2.state
    pp, rp, _, v = productcheck_prove(F, P, f1, f2, 3, transcript, label + b"/perm")
    return {"gate": pg, "perm": pp, "nonce": nonce}, v


def hyperplonk_verify(F, T, proof, v, transcript, label=b"hp"):
    """Verifier; T và v đóng vai oracle (truy vấn bằng mle.evaluate)."""
    from .mle import evaluate
    l = num_vars(T["a"])
    ok, r, final, eqv = zerocheck_verify(F, l, 3, proof["gate"], transcript, label + b"/gate")
    if not ok or eqv * _gate([evaluate(F, T[k], r) for k in GATE_KEYS]) != final:
        return False
    transcript.absorb(label + b"/nonce", int(proof["nonce"]))
    beta = _challenge(F, transcript, label + b"/beta")
    gamma = _challenge(F, transcript, label + b"/gamma")
    f1, f2 = _perm_fns(beta, gamma)
    ok, s = productcheck_verify(F, l, 3, proof["perm"], transcript, f1, f2,
                                lambda pt: evaluate(F, v, pt),
                                lambda pt: [evaluate(F, T[k], pt) for k in PERM_KEYS],
                                label + b"/perm")
    return bool(ok and s == 1)
