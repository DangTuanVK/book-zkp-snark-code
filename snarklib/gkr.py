"""Giao thức GKR cho mạch số học phân tầng (cổng cộng/nhân, fan-in 2).

Quy ước theo thaler2022proofs.txt:2876-2900 (Mục 4.6.4): tầng 0 là tầng ra,
tầng d là tầng vào; tầng i có S_i = 2^{k_i} cổng, W_i: {0,1}^{k_i} -> F là giá
trị cổng. Vị từ nối dây add_i, mult_i: {0,1}^{k_i + 2k_{i+1}} -> {0,1},
add_i(a, b, c) = 1 khi cổng a ở tầng i là cổng cộng với hai đầu vào b, c.
Bảng vị từ đánh chỉ số big-endian theo (a, b, c) nối liền (a là các bit cao).

Rút gọn tầng (Lemma 4.7, thaler2022proofs.txt:3002-3005, phương trình (4.17)):
    W_i~(z) = Σ_{b,c ∈ {0,1}^{k_{i+1}}} add_i~(z,b,c)·(W_{i+1}~(b) + W_{i+1}~(c))
                                       + mult_i~(z,b,c)·W_{i+1}~(b)·W_{i+1}~(c).
Khẳng định W_i~(r_i) = v được rút về hai khẳng định W_{i+1}~(b*) = z1,
W_{i+1}~(c*) = z2 bằng sum-check 2k_{i+1} biến, bậc 2 mỗi biến (dùng
snarklib.sumcheck với 4 bảng ADD, MUL, WB, WC). Gộp hai khẳng định bằng
đường thẳng ℓ(t) = b* + t(c* - b*) (thaler2022proofs.txt:3090-3097): Prover gửi
q(t) = W_{i+1}~(ℓ(t)) bậc ≤ k_{i+1}; Verifier kiểm q(0) = z1, q(1) = z2, rút
r*, khẳng định mới W_{i+1}~(ℓ(r*)) = q(r*). Tầng vào: Verifier tự tính W_d~.

Prover ở đây dựng bảng 2^{2k_{i+1}} ô (bậc hai theo số cổng); prover tuyến
tính (Libra, xie2019libra) tách sum-check thành hai pha — chỉ mô tả trong sách.
Verifier tự đánh giá add_i~, mult_i~ bằng bảng (O(S_i·S_{i+1}^2)), chấp nhận
được cho mạch nhỏ; mạch có cấu trúc đều đánh giá được trong O(k) (Mục 4.6.6).

API
    circ = LayeredCircuit(layers, n_inputs)   layers[i] = [("add"|"mul", in1, in2), ...]
    circ.evaluate(x) -> [W_0, W_1, ..., W_d]
    prove(circ, x, transcript)  -> proof
    verify(circ, x, outputs, proof, transcript) -> bool
"""

from .field import GF
from .mle import eq_table, evaluate, num_vars
from . import sumcheck as sc

F449 = GF(449)


def _pad_pow2(n):
    s = 1
    while s < n:
        s *= 2
    return s


class LayeredCircuit:
    """Mạch phân tầng. layers[0] là tầng ra; layers[-1] nối vào tầng vào (n_inputs)."""

    def __init__(self, layers, n_inputs, F=F449):
        self.F = F
        self.layers = [list(L) for L in layers]
        self.d = len(layers)                       # tầng d là tầng vào
        sizes = [_pad_pow2(len(L)) for L in layers] + [_pad_pow2(n_inputs)]
        self.sizes = sizes
        self.k = [num_vars([0] * s) for s in sizes]
        for i, L in enumerate(self.layers):
            for typ, a, b in L:
                assert typ in ("add", "mul") and a < sizes[i + 1] and b < sizes[i + 1]

    def evaluate(self, x):
        F = self.F
        W = [F(v) for v in x] + [F(0)] * (self.sizes[-1] - len(x))
        out = [W]
        for i in range(self.d - 1, -1, -1):
            nxt = []
            for typ, a, b in self.layers[i]:
                nxt.append(W[a] + W[b] if typ == "add" else W[a] * W[b])
            nxt += [F(0)] * (self.sizes[i] - len(nxt))
            W = nxt
            out.append(W)
        return out[::-1]                           # out[i] = W_i

    def predicate_table(self, i, typ):
        """Bảng add_i hoặc mult_i trên {0,1}^{k_i + 2k_{i+1}} (để Verifier đánh giá MLE)."""
        S, S1 = self.sizes[i], self.sizes[i + 1]
        T = [0] * (S * S1 * S1)
        for a, (t, b, c) in enumerate(self.layers[i]):
            if t == typ:
                T[(a * S1 + b) * S1 + c] = 1
        return T

    def predicate_eval(self, i, typ, z, bstar, cstar):
        return evaluate(self.F, self.predicate_table(i, typ), list(z) + list(bstar) + list(cstar))

    def wiring_tables(self, i, r):
        """ADD[(b,c)] = add_i~(r, b, c), MUL tương tự — tính thưa, O(S_i)."""
        F = self.F
        S1 = self.sizes[i + 1]
        er = eq_table(F, r)
        ADD = [F(0)] * (S1 * S1)
        MUL = [F(0)] * (S1 * S1)
        for a, (t, b, c) in enumerate(self.layers[i]):
            tgt = ADD if t == "add" else MUL
            tgt[b * S1 + c] = tgt[b * S1 + c] + er[a]
        return ADD, MUL


def _g(u):
    """add·(wb + wc) + mul·wb·wc  với u = (ADD, MUL, WB, WC)."""
    return u[0] * (u[2] + u[3]) + u[1] * u[2] * u[3]


def _line(F, bstar, cstar, t):
    t = F(t)
    return [b + t * (c - b) for b, c in zip(bstar, cstar)]


def prove(circ, x, transcript, label=b"gkr"):
    F = circ.F
    W = circ.evaluate(x)
    transcript.absorb(label + b"/out", [int(v) for v in W[0]])
    r = [sc._challenge(F, transcript, label + b"/r0_%d" % j) for j in range(circ.k[0])]
    v = evaluate(F, W[0], r)
    layers = []
    for i in range(circ.d):
        S1, k1 = circ.sizes[i + 1], circ.k[i + 1]
        ADD, MUL = circ.wiring_tables(i, r)
        WB = [W[i + 1][b] for b in range(S1) for _ in range(S1)]
        WC = [W[i + 1][c] for _ in range(S1) for c in range(S1)]
        lab = label + b"/L%d" % i
        proof, rr, fin = sc.prove(F, [ADD, MUL, WB, WC], _g, 2, transcript, claim=v, label=lab)
        bstar, cstar = rr[:k1], rr[k1:]
        z1, z2 = fin[2], fin[3]
        q = [evaluate(F, W[i + 1], _line(F, bstar, cstar, t)) for t in range(k1 + 1)]
        transcript.absorb(lab + b"/z", [int(z1), int(z2)])
        transcript.absorb(lab + b"/q", [int(e) for e in q])
        rstar = sc._challenge(F, transcript, lab + b"/rstar")
        r = _line(F, bstar, cstar, rstar)
        v = sc.interpolate_eval(F, q, rstar)
        layers.append({"sumcheck": proof, "z": (z1, z2), "q": q, "claim": proof["claim"]})
    return {"outputs": W[0], "layers": layers}


def verify(circ, x, outputs, proof, transcript, label=b"gkr", trace=None):
    """trace (list, tùy chọn) nhận (tầng, r_i, v_i) để in minh họa."""
    F = circ.F
    outs = [F(v) for v in outputs] + [F(0)] * (circ.sizes[0] - len(outputs))
    transcript.absorb(label + b"/out", [int(v) for v in outs])
    r = [sc._challenge(F, transcript, label + b"/r0_%d" % j) for j in range(circ.k[0])]
    v = evaluate(F, outs, r)
    if len(proof["layers"]) != circ.d:
        return False
    for i, L in enumerate(proof["layers"]):
        if trace is not None:
            trace.append((i, list(r), v))
        k1 = circ.k[i + 1]
        lab = label + b"/L%d" % i
        ok, rr, final = sc.verify(F, 2 * k1, 2, L["sumcheck"], transcript, lab, claim=v)
        if not ok:
            return False
        bstar, cstar = rr[:k1], rr[k1:]
        z1, z2 = F(L["z"][0]), F(L["z"][1])
        add_v = circ.predicate_eval(i, "add", r, bstar, cstar)
        mul_v = circ.predicate_eval(i, "mul", r, bstar, cstar)
        if add_v * (z1 + z2) + mul_v * z1 * z2 != final:
            return False
        q = [F(e) for e in L["q"]]
        if len(q) != k1 + 1 or q[0] != z1 or q[1] != z2:
            return False
        transcript.absorb(lab + b"/z", [int(z1), int(z2)])
        transcript.absorb(lab + b"/q", [int(e) for e in q])
        rstar = sc._challenge(F, transcript, lab + b"/rstar")
        r = _line(F, bstar, cstar, rstar)
        v = sc.interpolate_eval(F, q, rstar)
    xin = [F(e) for e in x] + [F(0)] * (circ.sizes[-1] - len(x))
    if trace is not None:
        trace.append((circ.d, list(r), v))
    return evaluate(F, xin, r) == v


def proof_size(circ):
    """Số phần tử trường trong proof (không kể đầu ra): mỗi tầng 2k·3 + 2 + (k+1)."""
    return sum(2 * circ.k[i + 1] * 3 + 2 + circ.k[i + 1] + 1 for i in range(circ.d))


def vi_du_mach():
    """Mạch 3 tầng của sách: vào (a, b, c, d) = (17, 23, 28, 30).

    Tầng 1: [a·b, c·d, a+b, c+d] = [391, 840 ≡ 391, 40, 58]
    Tầng 0: [a·b + c·d, (a+b)·(c+d)] = [782 ≡ 333, 2320 ≡ 75]   (mod 449)
    """
    layers = [
        [("add", 0, 1), ("mul", 2, 3)],                                   # tầng 0
        [("mul", 0, 1), ("mul", 2, 3), ("add", 0, 1), ("add", 2, 3)],     # tầng 1
    ]
    return LayeredCircuit(layers, 4), [17, 23, 28, 30]
