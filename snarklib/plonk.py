"""mini-PLONK đầy đủ (Chương 6): 5 vòng Prover và 12 bước Verifier.

Theo đúng thứ tự của bản ePrint 2019/953 "PlonK: Permutations over
Lagrange-bases for Oecumenical Noninteractive arguments of Knowledge"
(Gabizon, Williamson, Ciobotaru), bản ngày 17/07/2025, file
book-zkp-snark/tai-lieu-tham-khao/txt/gabizon2019plonk.txt:
    §8.1 đa thức của mạch ............................ dòng 1250-1278
    §8.3 Common preprocessed input ................... dòng 1295-1311
    Round 1 (đa thức dây có làm mù b1..b6) ........... dòng 1321-1341
    Round 2 (β, γ và z(X), làm mù b7..b9) ............ dòng 1345-1363
    Round 3 (α và t(X), tách t_lo/t_mid/t_hi, b10,b11) dòng 1369-1393
    Round 4 (ζ và 6 giá trị mở) ...................... dòng 1397-1409
    Round 5 (v, r(X), W_ζ, W_ζω, sau đó u) ........... dòng 1413-1462
    Verifier 12 bước ................................. dòng 1466-1529
    Chi phí: Bảng 1, 2 (dòng 140-181): 9 G1 + 6 F, 2 pairing, 18 lũy thừa G1.

CHỖ LỆCH SO VỚI BÀI (ghi rõ để văn bản sách theo):
  (1) Chỉ số từ 0: hàng j ứng với ω^j, j = 0..n-1 (bài dùng i = 1..n và ω^i).
      Vì vậy L_1 của bài (ràng buộc khởi đầu z = 1 và PI) thành L_0 ở đây,
      L_0(ζ) = (ζ^n - 1)/(n(ζ - 1)) thay cho ω(ζ^n - 1)/(n(ζ - ω)) ở bước 6.
  (2) PI(X) = Σ -x_i L_i(X) theo §5 của bài (dòng 1164) và quy ước §0.4 của
      sách (q_L = 1, PI_i = -x_i). Bước 7 của Verifier trong txt viết
      PI(z) = Σ w_i L_i(z) (không có dấu trừ) -- mâu thuẫn dấu trong chính bài;
      ta dùng dấu trừ cho khớp phương trình cổng.
  (3) Thách thức của bài ký hiệu z (điểm đánh giá) trùng tên với z(X); ở đây
      gọi là ζ. Thách thức β = H(transcript, 0), γ = H(transcript, 1) được cài
      bằng hai lần challenge liên tiếp với nhãn khác nhau của snarklib.Transcript.
  (4) Vòng 3 tính t(X) bằng coset FFT. Bài (và các cài đặt thực tế) đánh giá
      trên một coset cỡ 4n; F_449 chỉ có căn đơn vị bậc ≤ 64 nên với n = 32
      không có miền cỡ 128. Ta dùng 4 coset s_j·H cỡ n (mỗi coset một NTT cỡ n)
      và ghép lại bằng định lý phần dư Trung Hoa theo X^n - s_j^n (xem
      interpolate_from_cosets). Kết quả trùng với phép chia đa thức trực tiếp.
  (5) Mạch có 20 hàng dùng được, đệm hàng 0 tới n = 32; ô không dùng nhận giá
      trị 0 và không bị ràng buộc sao chép (σ là đồng nhất tại đó).

Cách dùng:
    pk, vk = setup(circuit, tau)                      # circuit: PlonkishCircuit
    proof, log = prove(pk, a, b, c, x, rng)           # a, b, c: cột dây (≤ n hàng)
    ok = verify(vk, x, proof)
Tham số hash_public=False ở prove/verify tái hiện lỗi "Frozen Heart" (không
băm đầu vào công khai vào transcript); xem forge_frozen_heart.
"""

import random
from collections import namedtuple

from .field import GF
from .poly import Poly, vanishing_poly, lagrange_interpolate
from .fft import intt, coset_ntt, coset_intt
from .transcript import Transcript
from .curve import R_TOY
from . import kzg
from . import opcount as oc

R = R_TOY
F = GF(R)
K1, K2 = 2, 3                 # coset H, k1·H, k2·H (hang_so.txt)
N_COSETS = 4                  # số coset cỡ n dùng ở vòng 3 (4n > 3n + 5)

Proof = namedtuple("Proof", [
    "a", "b", "c", "z", "t_lo", "t_mid", "t_hi", "W_zeta", "W_zeta_omega",   # 9 G1
    "a_bar", "b_bar", "c_bar", "s1_bar", "s2_bar", "z_omega_bar"])          # 6 F
G1_FIELDS = Proof._fields[:9]
FR_FIELDS = Proof._fields[9:]


class VerifyingKey:
    """Verifier preprocessed input (dòng 1470-1473): 8 cam kết + [1]_1, [1]_2, [x]_2."""

    def __init__(self, n, omega, ell, k1, k2, comms, G1, g2):
        self.n, self.omega, self.ell, self.k1, self.k2 = n, omega, ell, k1, k2
        self.comms = comms        # dict: qM qL qR qO qC s1 s2 s3 -> điểm G1
        self.G1, self.g2 = G1, g2  # g2 = [[1]_2, [x]_2]


class ProvingKey:
    """Common preprocessed input (dòng 1297-1311) + bảng đánh giá trên coset."""

    def __init__(self, **kw):
        self.__dict__.update(kw)


# ---------------------------------------------------------------------------
# Đếm NTT (opcount chỉ có sẵn phép toán nhóm; ta ghi thêm khóa "ntt")
# ---------------------------------------------------------------------------

def _count_ntt(n):
    oc._bump("ntt", 1)
    oc._bump("ntt_size", n)


def _intt(vals):
    _count_ntt(len(vals))
    return intt(F, vals)


def _ntt_coset(coeffs, s):
    _count_ntt(len(coeffs))
    return coset_ntt(F, coeffs, s)


def eval_on_coset(f, s, n):
    """[f(s·ω^i)]_{i<n} bằng MỘT NTT cỡ n: rút gọn f mod (X^n - s^n) trước.

    Trên coset sH ta có X^n = s^n nên hệ số bậc k + m·n cộng vào bậc k với
    trọng số (s^n)^m; bậc của f có thể lớn hơn n (a, b, c, z có bậc n+1, n+2).
    """
    s = F(s)
    sn, folded, pw = s ** n, [F(0)] * n, F(1)
    for m in range(0, len(f.coeffs), n):
        for k, cf in enumerate(f.coeffs[m:m + n]):
            folded[k] = folded[k] + cf * pw
        pw = pw * sn
    return _ntt_coset(folded, s)


def interpolate_from_cosets(vals, shifts, n):
    """Khôi phục đa thức bậc < len(shifts)·n từ giá trị trên các coset s_j·H.

    coset_intt trên s_j·H cho R_j = t mod (X^n - c_j), c_j = s_j^n. Viết
    t = Σ_m X^{mn} T_m với deg T_m < n thì R_j = Σ_m c_j^m T_m: với mỗi vị trí
    hệ số k, (T_m[k])_m là hệ số của đa thức theo biến Y nội suy qua các điểm
    (c_j, R_j[k]) -- phép nội suy Lagrange cỡ len(shifts).
    """
    Rs = []
    for v, s in zip(vals, shifts):
        _count_ntt(n)
        Rs.append(coset_intt(F, v, s))
    cs = [F(s) ** n for s in shifts]
    out = [F(0)] * (n * len(shifts))
    for k in range(n):
        P = lagrange_interpolate(F, cs, [Rj[k] for Rj in Rs])
        for m, cf in enumerate(P.coeffs):
            out[m * n + k] = cf
    return Poly(F, out)


def coset_shifts(n, count=N_COSETS):
    """count phần tử s ≥ 5 nhỏ nhất có s^n ≠ 1 và các s^n đôi một khác nhau.

    Bắt đầu từ 5 để coset không trùng k1·H, k2·H (cho dễ đọc; về toán không bắt buộc).
    """
    out, seen, s = [], set(), 5
    while len(out) < count:
        c = int(F(s) ** n)
        if c != 1 and c not in seen and c not in (int(F(K1) ** n), int(F(K2) ** n)):
            out.append(F(s))
            seen.add(c)
        s += 1
    return out


# ---------------------------------------------------------------------------
# Tiền xử lý: selector, hoán vị σ*, SRS
# ---------------------------------------------------------------------------

def domain_size(rows):
    n = 1
    while n < rows:
        n *= 2
    return n


def build_sigma(copies, n, omega, k1=K1, k2=K2):
    """σ* (3n giá trị) từ danh sách ràng buộc sao chép ((cột, hàng), (cột, hàng)).

    Ô (cột c, hàng j) đánh số c·n + j và mang nhãn đồng nhất id = k_c·ω^j
    (k_a = 1, k_b = k1, k_c = k2). Các ô bằng nhau gom thành lớp (union-find);
    trong mỗi lớp, σ chuyển mỗi ô sang ô kế tiếp theo vòng. σ*(ô) = id(σ(ô)).
    Trả về (sigma, ids, cycles): sigma là hoán vị chỉ số trên [3n].
    """
    parent = list(range(3 * n))

    def find(u):
        while parent[u] != u:
            parent[u] = parent[parent[u]]
            u = parent[u]
        return u

    col = {"a": 0, "b": 1, "c": 2}
    for (c1, r1), (c2, r2) in copies:
        u, v = find(col[c1] * n + r1), find(col[c2] * n + r2)
        if u != v:
            parent[max(u, v)] = min(u, v)
    classes = {}
    for u in range(3 * n):
        classes.setdefault(find(u), []).append(u)
    sigma = list(range(3 * n))
    cycles = []
    for cls in classes.values():
        if len(cls) > 1:
            cycles.append(cls)
        for i, u in enumerate(cls):
            sigma[u] = cls[(i + 1) % len(cls)]
    ks = [F(1), F(k1), F(k2)]
    ids = [ks[u // n] * omega ** (u % n) for u in range(3 * n)]
    return sigma, ids, cycles


def setup(circuit, tau, n=None, k1=K1, k2=K2):
    """Sinh (pk, vk) cho mạch Plonkish; tau là cửa sập SRS (chỉ để minh họa).

    SRS: [x^i]_1, i = 0..n+5 (dòng 1298: "x^{n+5}·[1]_1") và [1]_2, [x]_2.
    Thời gian tiền xử lý không đếm vào opcount của Prover/Verifier.
    """
    rows = circuit.n
    n = n or domain_size(rows)
    omega = F.root_of_unity(n)
    for k in (k1, k2, F(k1) / F(k2)):                 # ba coset phân biệt
        assert F(k) ** n != 1, "k1, k2 không cho ba coset phân biệt"
    if F(tau) ** n == 1 or F(tau) == 0:
        raise ValueError("tau thuộc H (Z_H(tau) = 0): làm mù vô hiệu, cam kết suy biến")
    pad = lambda v: list(v) + [F(0)] * (n - len(v))
    with oc.scope():                                   # bỏ đếm phần tiền xử lý
        sel_vals = {s: pad(circuit.sel[s]) for s in ("qM", "qL", "qR", "qO", "qC")}
        sigma, ids, cycles = build_sigma(circuit.copies, n, omega, k1, k2)
        sstar = [ids[sigma[u]] for u in range(3 * n)]
        polys = {s: Poly(F, intt(F, v)) for s, v in sel_vals.items()}
        for j in range(3):
            polys["s%d" % (j + 1)] = Poly(F, intt(F, sstar[j * n:(j + 1) * n]))
        srs = kzg.setup(n + 5, tau)
        comms = {k: kzg.commit(srs, p) for k, p in polys.items()}
        shifts = coset_shifts(n)
        L0 = Poly(F, intt(F, [F(1)] + [F(0)] * (n - 1)))
        polys["L0"] = L0
        coset = {k: [eval_on_coset(p, s, n) for s in shifts] for k, p in polys.items()}
    vk = VerifyingKey(n, omega, circuit.n_public, F(k1), F(k2), comms, srs.G1, srs.g2)
    pk = ProvingKey(n=n, omega=omega, ell=circuit.n_public, k1=F(k1), k2=F(k2),
                    rows=rows, sel_vals=sel_vals, sigma=sigma, ids=ids, sstar=sstar,
                    cycles=cycles, polys=polys, srs=srs, vk=vk, shifts=shifts,
                    coset=coset, circuit=circuit)
    return pk, vk


def new_transcript(vk, x, label=b"snarklib/plonk", hash_public=True):
    """Transcript = common preprocessed input ‖ public input (§8.3, dòng 1283-1287).

    hash_public=False bỏ x khỏi transcript: đây chính là lỗi Fiat-Shamir mà bài
    cảm ơn Sean Bowe đã phát hiện (dòng 1528-1529) -- dùng cho demo giả mạo.
    """
    tr = Transcript(label, R)
    tr.absorb(b"n", vk.n)
    tr.absorb(b"k", [vk.k1, vk.k2])
    tr.absorb(b"preprocessed", [vk.comms[k] for k in
                                ("qM", "qL", "qR", "qO", "qC", "s1", "s2", "s3")])
    if hash_public:
        tr.absorb(b"public", [F(v) for v in x])
    return tr


def pi_poly(n, x):
    """PI(X) = Σ_{i<ℓ} -x_i L_i(X) (chỉ số từ 0)."""
    return Poly(F, intt(F, [-F(v) for v in x] + [F(0)] * (n - len(x))))


def lagrange0_at(n, zeta):
    """L_0(ζ) = (ζ^n - 1) / (n(ζ - 1)); tại ζ = 1 công thức có dạng 0/0 và L_0(1) = 1."""
    if zeta == F(1):
        return F(1)
    return (zeta ** n - 1) / (F(n) * (zeta - 1))


def pi_at(n, omega, x, zeta):
    """PI(ζ) = Σ -x_i L_i(ζ), L_i(ζ) = ω^i(ζ^n - 1)/(n(ζ - ω^i))."""
    zh = zeta ** n - 1
    return sum((-F(v) * omega ** i * zh / (F(n) * (zeta - omega ** i))
                for i, v in enumerate(x)), F(0))


# ---------------------------------------------------------------------------
# Prover: 5 vòng (mỗi vòng một phương thức; lớp con có thể gắn thêm cột/cổng)
# ---------------------------------------------------------------------------

class Prover:
    """Prover PLONK. log ghi (vòng, nhãn, giá trị) theo đúng thứ tự transcript."""

    def __init__(self, pk, a, b, c, x, rng=None, label=b"snarklib/plonk",
                 hash_public=True, check=True):
        self.pk, self.n, self.check = pk, pk.n, check
        pad = lambda v: [F(t) for t in v] + [F(0)] * (pk.n - len(v))
        self.w = [pad(a), pad(b), pad(c)]
        self.x = [F(v) for v in x]
        self.rng = rng or random.Random()
        self.tr = new_transcript(pk.vk, self.x, label, hash_public)
        self.log = []
        self.ZH = vanishing_poly(F, pk.n)

    def _log(self, rnd, name, val):
        self.log.append((rnd, name, val))

    def _commit(self, f):
        return kzg.commit(self.pk.srs, f)

    def _blinded(self, vals, blinds):
        """Σ w_i L_i(X) + (b_1 X^{k-1} + ... + b_k) Z_H(X) (dòng 1325-1336, 1353)."""
        return Poly(F, _intt(vals)) + Poly(F, list(reversed(blinds))) * self.ZH

    # Round 1 ---------------------------------------------------------------
    def round1(self):
        self.b = [None] + [F(self.rng.randrange(R)) for _ in range(11)]   # b1..b11
        b = self.b
        self.a_p = self._blinded(self.w[0], [b[1], b[2]])
        self.b_p = self._blinded(self.w[1], [b[3], b[4]])
        self.c_p = self._blinded(self.w[2], [b[5], b[6]])
        self.cm_abc = [self._commit(p) for p in (self.a_p, self.b_p, self.c_p)]
        self.tr.absorb(b"round1", self.cm_abc)
        for nm, cm in zip("abc", self.cm_abc):
            self._log(1, "[%s]_1" % nm, cm)

    # Round 2 ---------------------------------------------------------------
    def round2(self):
        pk, n = self.pk, self.n
        self.beta = self.tr.challenge(b"beta")     # H(transcript, 0)
        self.gamma = self.tr.challenge(b"gamma")   # H(transcript, 1)
        self._log(2, "beta", self.beta)
        self._log(2, "gamma", self.gamma)
        be, ga, ids, ss, w = self.beta, self.gamma, pk.ids, pk.sstar, self.w
        zv = [F(1)]
        for i in range(n - 1):
            num = den = F(1)
            for j in range(3):
                num = num * (w[j][i] + be * ids[j * n + i] + ga)
                den = den * (w[j][i] + be * ss[j * n + i] + ga)
            if den == 0:
                raise ZeroDivisionError("mẫu số z bằng 0: Prover dừng (chú thích 11 của bài)")
            zv.append(zv[-1] * num / den)
        self.z_vals = zv
        self.z_p = self._blinded(zv, [self.b[7], self.b[8], self.b[9]])
        self.cm_z = self._commit(self.z_p)
        self.tr.absorb(b"round2", self.cm_z)
        self._log(2, "[z]_1", self.cm_z)

    # Round 3 ---------------------------------------------------------------
    def numerator_on_coset(self, j, ev):
        """Tử số của t(X) trên coset thứ j; ev: dict tên -> list giá trị.

        Lớp con (ví dụ LogUp ở Ch.7) gọi super() rồi cộng thêm các hạng α^3, α^4...
        """
        be, ga, al = self.beta, self.gamma, self.alpha
        k1, k2 = self.pk.k1, self.pk.k2
        out = []
        for i in range(self.n):
            A, B, C, Z, Zw, X = (ev[k][i] for k in ("a", "b", "c", "z", "zw", "X"))
            gate = (A * B * ev["qM"][i] + A * ev["qL"][i] + B * ev["qR"][i]
                    + C * ev["qO"][i] + ev["PI"][i] + ev["qC"][i])
            perm1 = (A + be * X + ga) * (B + be * k1 * X + ga) * (C + be * k2 * X + ga) * Z
            perm2 = ((A + be * ev["s1"][i] + ga) * (B + be * ev["s2"][i] + ga)
                     * (C + be * ev["s3"][i] + ga) * Zw)
            out.append(gate + al * (perm1 - perm2) + al * al * (Z - 1) * ev["L0"][i])
        return out

    def quotient_coset_fft(self):
        """t(X) = tử số / Z_H bằng coset FFT trên N_COSETS coset s_j·H (vòng 3)."""
        pk, n = self.pk, self.n
        PI = pi_poly(n, self.x)
        tvals = []
        for j, s in enumerate(pk.shifts):
            ev = {k: pk.coset[k][j] for k in ("qM", "qL", "qR", "qO", "qC",
                                               "s1", "s2", "s3", "L0")}
            ev["a"], ev["b"], ev["c"] = (eval_on_coset(p, s, n)
                                         for p in (self.a_p, self.b_p, self.c_p))
            ev["z"] = eval_on_coset(self.z_p, s, n)
            ev["zw"] = ev["z"][1:] + ev["z"][:1]          # z(ω·sω^i) = z(sω^{i+1})
            ev["PI"] = eval_on_coset(PI, s, n)
            ev["X"] = [s * pk.omega ** i for i in range(n)]
            ev.update(self.extra_coset_evals(j, s))
            zh_inv = (s ** n - 1).inverse()               # Z_H hằng trên coset
            tvals.append([v * zh_inv for v in self.numerator_on_coset(j, ev)])
        return interpolate_from_cosets(tvals, pk.shifts, n)

    def extra_coset_evals(self, j, s):
        """Móc mở rộng: lớp con trả thêm giá trị trên coset (cột mới)."""
        return {}

    def quotient_naive(self):
        """Cùng t(X) bằng phép chia đa thức (đối chiếu; dùng khi check=False)."""
        P, n = self.pk.polys, self.n
        be, ga, al, k1, k2 = self.beta, self.gamma, self.alpha, self.pk.k1, self.pk.k2
        X = Poly.X(F)
        a, b, c, z = self.a_p, self.b_p, self.c_p, self.z_p
        zw = z.scale(self.pk.omega)
        num = (a * b * P["qM"] + a * P["qL"] + b * P["qR"] + c * P["qO"]
               + pi_poly(n, self.x) + P["qC"])
        num = num + ((a + X * be + ga) * (b + X * (be * k1) + ga)
                     * (c + X * (be * k2) + ga) * z) * al
        num = num - ((a + P["s1"] * be + ga) * (b + P["s2"] * be + ga)
                     * (c + P["s3"] * be + ga) * zw) * al
        num = num + ((z - 1) * P["L0"]) * (al * al)
        return num.divmod(self.ZH)

    def round3(self, t_shift=0):
        """t_shift ≠ 0 chỉ dùng cho demo giả mạo (cộng hằng số vào t_lo)."""
        n = self.n
        self.alpha = self.tr.challenge(b"alpha")
        self._log(3, "alpha", self.alpha)
        if self.check:
            t = self.quotient_coset_fft()
            if t.degree() > 3 * n + 5:
                raise ValueError("witness không thỏa mạch: tử số không chia hết cho Z_H")
        else:
            t, _ = self.quotient_naive()      # Prover gian: bỏ phần dư
        self.t_p = t
        co = t.coeffs + [F(0)] * (3 * n + 6 - len(t.coeffs))
        lo, mid, hi = Poly(F, co[:n]), Poly(F, co[n:2 * n]), Poly(F, co[2 * n:])
        b10, b11 = self.b[10], self.b[11]
        Xn = Poly(F, [0] * n + [1])
        self.t_lo = lo + Xn * b10 + t_shift
        self.t_mid = mid - b10 + Xn * b11
        self.t_hi = hi - b11
        self._commit_t()

    def _commit_t(self):
        self.cm_t = [self._commit(p) for p in (self.t_lo, self.t_mid, self.t_hi)]
        self.tr_before_zeta = _copy_tr(self.tr)
        self.tr.absorb(b"round3", self.cm_t)
        for nm, cm in zip(("lo", "mid", "hi"), self.cm_t):
            self._log(3, "[t_%s]_1" % nm, cm)

    # Round 4 ---------------------------------------------------------------
    def round4(self):
        P, z = self.pk.polys, self.tr.challenge(b"zeta")
        self.zeta = z
        self._log(4, "zeta", z)
        self.evals = [self.a_p(z), self.b_p(z), self.c_p(z), P["s1"](z), P["s2"](z),
                      self.z_p(z * self.pk.omega)]
        self.tr.absorb(b"round4", self.evals)
        for nm, v in zip(FR_FIELDS, self.evals):
            self._log(4, nm, v)

    # Round 5 ---------------------------------------------------------------
    def linearization(self):
        """r(X) (dòng 1419-1426); r(ζ) = 0 khi mọi ràng buộc thỏa."""
        P, n, z = self.pk.polys, self.n, self.zeta
        be, ga, al, k1, k2 = self.beta, self.gamma, self.alpha, self.pk.k1, self.pk.k2
        a, b, c, s1, s2, zw = self.evals
        r = (P["qM"] * (a * b) + P["qL"] * a + P["qR"] * b + P["qO"] * c
             + pi_at(n, self.pk.omega, self.x, z) + P["qC"])
        r = r + self.z_p * (al * (a + be * z + ga) * (b + be * k1 * z + ga)
                            * (c + be * k2 * z + ga))
        r = r - (P["s3"] * be + (c + ga)) * (al * (a + be * s1 + ga) * (b + be * s2 + ga) * zw)
        r = r + (self.z_p - 1) * (al * al * lagrange0_at(n, z))
        r = r - (self.t_lo + self.t_mid * z ** n + self.t_hi * z ** (2 * n)) * (z ** n - 1)
        return r

    def zeta_openings(self):
        """Danh sách (đa thức, giá trị) gộp vào W_ζ theo lũy thừa v^0, v^1, ...

        Lớp con gắn thêm cột thì nối thêm cặp vào cuối (và Verifier làm tương ứng).
        """
        P, (a, b, c, s1, s2, _) = self.pk.polys, self.evals
        return [(self.linearization(), F(0)), (self.a_p, a), (self.b_p, b),
                (self.c_p, c), (P["s1"], s1), (P["s2"], s2)]

    def round5(self):
        z, om = self.zeta, self.pk.omega
        self.v = self.tr.challenge(b"v")
        self._log(5, "v", self.v)
        acc, vp = Poly.zero(F), F(1)
        for f, y in self.zeta_openings():
            acc = acc + (f - y) * vp
            vp = vp * self.v
        Wz, rem = acc.divmod(Poly(F, [-z, 1]))
        Wzw, rem2 = (self.z_p - self.evals[5]).divmod(Poly(F, [-z * om, 1]))
        if self.check:
            assert rem.is_zero() and rem2.is_zero(), "r(ζ) ≠ 0"
        self.cm_W = [self._commit(Wz), self._commit(Wzw)]
        self.tr.absorb(b"round5", self.cm_W)
        self._log(5, "[W_zeta]_1", self.cm_W[0])
        self._log(5, "[W_zeta_omega]_1", self.cm_W[1])
        self.u = self.tr.challenge(b"u")            # multipoint evaluation challenge
        self._log(5, "u", self.u)

    def proof(self):
        return Proof(*(self.cm_abc + [self.cm_z] + self.cm_t + self.cm_W + self.evals))

    def run(self, max_aborts=50):
        """Chạy 5 vòng. Nếu mẫu số của z bằng 0 (xác suất ≈ 3(n-1)/r, khoảng 19%
        trên F_449 với n = 32, không đáng kể ở trường 256 bit), Prover rút lại
        hệ số làm mù b1..b11 và chạy lại từ vòng 1; số lần dừng ghi ở self.aborts."""
        self.aborts = 0
        tr0 = _copy_tr(self.tr)
        while True:
            self.round1()
            try:
                self.round2()
                break
            except ZeroDivisionError:
                self.aborts += 1
                if self.aborts >= max_aborts:
                    raise
                self.tr, self.log = _copy_tr(tr0), []
        self.round3()
        self.round4()
        self.round5()
        return self.proof()


def _copy_tr(tr):
    t = Transcript.__new__(Transcript)
    t.F, t.state = tr.F, tr.state
    return t


def prove(pk, a, b, c, x, rng=None, label=b"snarklib/plonk", hash_public=True,
          check=True):
    """Trả về (proof, log). check=False: Prover gian (witness sai vẫn ra proof)."""
    p = Prover(pk, a, b, c, x, rng, label, hash_public, check)
    proof = p.run()
    return proof, p.log


# ---------------------------------------------------------------------------
# Verifier: 12 bước đúng thứ tự dòng 1477-1529
# ---------------------------------------------------------------------------

def _in_group(P):
    return P.on_curve() and (R * P).is_infinity()   # kiểm nhóm con (không đếm opcount)


def challenges(vk, x, proof, label=b"snarklib/plonk", hash_public=True):
    """Bước 4: β, γ, α, ζ, v, u từ transcript (giống hệt Prover)."""
    tr = new_transcript(vk, x, label, hash_public)
    tr.absorb(b"round1", [proof.a, proof.b, proof.c])
    be, ga = tr.challenge(b"beta"), tr.challenge(b"gamma")
    tr.absorb(b"round2", proof.z)
    al = tr.challenge(b"alpha")
    tr.absorb(b"round3", [proof.t_lo, proof.t_mid, proof.t_hi])
    ze = tr.challenge(b"zeta")
    tr.absorb(b"round4", [getattr(proof, k) for k in FR_FIELDS])
    v = tr.challenge(b"v")
    tr.absorb(b"round5", [proof.W_zeta, proof.W_zeta_omega])
    u = tr.challenge(b"u")
    return be, ga, al, ze, v, u


def verify(vk, x, proof, label=b"snarklib/plonk", hash_public=True, trace=None):
    """Verifier PLONK; trace (dict) nhận các giá trị trung gian để in sách."""
    n, om, k1, k2, cm = vk.n, vk.omega, vk.k1, vk.k2, vk.comms
    # 1-3. Kiểm tính hợp lệ của 9 phần tử G1, 6 phần tử F và đầu vào công khai
    if not all(_in_group(getattr(proof, k)) for k in G1_FIELDS):
        return False
    a, b, c, s1, s2, zw = (F(int(getattr(proof, k))) for k in FR_FIELDS)
    if len(x) != vk.ell:
        return False
    x = [F(v) for v in x]
    # 4. Thách thức
    be, ga, al, ze, v, u = challenges(vk, x, proof, label, hash_public)
    # 5. Z_H(ζ) = ζ^n - 1      6. L_0(ζ)      7. PI(ζ)
    zh = ze ** n - 1
    l0 = lagrange0_at(n, ze)
    pi = pi_at(n, om, x, ze)
    # 8. Hằng số r_0 của r(X)
    r0 = pi - l0 * al * al - al * (a + be * s1 + ga) * (b + be * s2 + ga) * (c + ga) * zw
    # 9. [D]_1 = [r']_1 + u·[z]_1 (hệ số của từng cam kết)
    sc_z = (a + be * ze + ga) * (b + be * k1 * ze + ga) * (c + be * k2 * ze + ga) * al \
        + l0 * al * al + u
    sc_s3 = -(a + be * s1 + ga) * (b + be * s2 + ga) * al * be * zw
    D = [(cm["qM"], a * b), (cm["qL"], a), (cm["qR"], b), (cm["qO"], c), (cm["qC"], F(1)),
         (proof.z, sc_z), (cm["s3"], sc_s3), (proof.t_lo, -zh),
         (proof.t_mid, -zh * ze ** n), (proof.t_hi, -zh * ze ** (2 * n))]
    # 10. [F]_1 = [D]_1 + v[a] + v^2[b] + v^3[c] + v^4[s1] + v^5[s2]
    Fl = D + [(proof.a, v), (proof.b, v ** 2), (proof.c, v ** 3),
              (cm["s1"], v ** 4), (cm["s2"], v ** 5)]
    # 11. [E]_1 = (-r_0 + vā + v^2 b̄ + v^3 c̄ + v^4 s̄1 + v^5 s̄2 + u z̄ω)·[1]_1
    e = -r0 + v * a + v ** 2 * b + v ** 3 * c + v ** 4 * s1 + v ** 5 * s2 + u * zw
    # 12. e([W_ζ] + u[W_ζω], [x]_2) = e(ζ[W_ζ] + uζω[W_ζω] + [F] - [E], [1]_2)
    #     Vế phải gộp thành MỘT MSM cỡ 18 (17 vô hướng khác 1), vế trái MSM cỡ 2.
    pts = [P for P, _ in Fl] + [vk.G1, proof.W_zeta, proof.W_zeta_omega]
    scs = [s for _, s in Fl] + [-e, ze, u * ze * om]
    rhs = oc.msm(pts, scs)
    lhs = oc.msm([proof.W_zeta, proof.W_zeta_omega], [F(1), u])
    ok = oc.pairing_check([(lhs, vk.g2[1]), (-rhs, vk.g2[0])])
    if trace is not None:
        trace.update(dict(beta=be, gamma=ga, alpha=al, zeta=ze, v=v, u=u, ZH=zh,
                          L0=l0, PI=pi, r0=r0, E=e, ok=ok))
    return ok


# ---------------------------------------------------------------------------
# Demo Frozen Heart: giả mạo khi transcript KHÔNG băm đầu vào công khai
# ---------------------------------------------------------------------------

def forge_frozen_heart(pk, a, b, c, x_honest, target=None, rng=None,
                       label=b"snarklib/plonk", hash_public=False, max_tries=20000):
    """Giả mạo proof PLONK cho đầu vào công khai sai (ℓ = 1).

    CƠ CHẾ. Verifier chỉ dùng x qua hằng số PI(ζ) = -x·L_0(ζ) trong r_0 (bước
    7-8). Nếu x không nằm trong transcript, mọi thách thức β, γ, α, ζ, v, u chỉ
    phụ thuộc các cam kết, nên Prover có thể chọn x SAU KHI biết ζ:
      1. Chạy vòng 1-2 trung thực với witness cho x_honest (391 = 17·23).
      2. Vòng 3: t(X) = t_honest(X) + Q(X), Q = q0 + q1·X cộng vào t_lo. Đẳng
         thức tại ζ giờ lệch Z_H(ζ)·Q(ζ), tức cần PI'(ζ) = PI_391(ζ) + Z_H(ζ)Q(ζ).
      3. Nhận ζ; vì Z_H(ζ)/L_0(ζ) = n(ζ - 1) nên đầu vào khớp là
             x* = x_honest - n(ζ - 1)·Q(ζ).
      4. Vòng 4-5 tính với PI'(ζ) của x* (để r(ζ) = 0, W_ζ chia hết). Proof
         verify với x* dù Prover không có witness cho x*.
    Với target = None: một lần thử (q0 = 1) cho proof hợp lệ của MỘT x* (không
    chọn trước) -- xác suất thành công 1, đúng cả ở trường 256 bit.
    Với target cho trước (vd 389): x* = target là sự kiện xác suất ≈ 1/r mỗi lần
    thử, nên phải "mài" (q0, q1) = (1,0), (2,0), ...; khả thi chỉ vì r = 449 nhỏ. Lưu ý: trên
    trường toy, việc mài này cũng phá được bản CÓ băm x (hash_public=True: x cố
    định trước, mỗi q0 cho ζ mới) -- đó là sai số soundness ~1/449 của trường
    toy, không phải lỗi Frozen Heart. Ở kích thước thật chỉ bản không băm bị phá.
    Trả về (proof, x_star, so_lan_thu) hoặc (None, None, max_tries).
    """
    n = pk.n
    p = Prover(pk, a, b, c, x_honest if not hash_public else [target], rng, label,
               hash_public, check=False)
    if hash_public:
        p.x = [F(v) for v in x_honest]     # t trung thực theo 391, transcript băm target
    tr0 = _copy_tr(p.tr)
    while True:
        p.round1()
        try:
            p.round2()
            break
        except ZeroDivisionError:
            p.tr, p.log = _copy_tr(tr0), []
    p.round3(t_shift=0)
    base_lo, (cm0, cm_mid, cm_hi) = p.t_lo, p.cm_t
    base_tr = _copy_tr(p.tr_before_zeta)
    G, Gt = pk.srs.g1[0], pk.srs.g1[1]
    for tries in range(1, max_tries + 1):
        # Q(X) = q0 + q1·X; q0 chạy 1..r-1, hết vòng thì tăng q1 (tránh lặp chu kỳ r)
        q1, q0 = F((tries - 1) // (R - 1)), F((tries - 1) % (R - 1) + 1)
        cm_lo = cm0 + int(q0) * G + int(q1) * Gt      # [t_lo + Q]_1 (không đếm opcount)
        t = _copy_tr(base_tr)
        t.absorb(b"round3", [cm_lo, cm_mid, cm_hi])
        zeta = t.challenge(b"zeta")
        x_star = F(x_honest[0]) - F(n) * (zeta - 1) * (q0 + q1 * zeta)
        if target is None or x_star == F(target):
            p.tr = _copy_tr(base_tr)
            p.log = [e for e in p.log if e[0] < 3 or e[1] == "alpha"]
            p.t_lo = base_lo + Poly(F, [q0, q1])
            p._commit_t()
            p.x = [x_star]
            p.round4()
            p.round5()
            p.forge_q = (int(q0), int(q1))
            return p.proof(), int(x_star), tries
    return None, None, max_tries
