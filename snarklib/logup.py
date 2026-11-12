"""LogUp: lookup bằng đạo hàm logarit (Chương 7, §7.2a và §7.3).

Nguồn: U. Haböck, "Multivariate lookups based on logarithmic derivatives",
ePrint 2022/1530 -- CHƯA có bản txt trong tai-lieu-tham-khao/txt/, mọi chi
tiết dưới đây cài theo hiểu biết tốt nhất và đang chờ đối chiếu (mục
pending_source_checks của báo cáo đợt 2b, ket_qua/bao_cao_dot/).

1. Đồng nhất thức (eq:ch7-logup). Lấy đạo hàm logarit hai vế của
   Π_i (X + f_i) = Π_j (X + t_j)^{m_j} (snarklib.multiset, (2)):

        Σ_i 1/(X + f_i)  =  Σ_j m_j/(X + t_j)                               (3)

   Bội m_j giờ là HỆ SỐ trong F_p chứ không còn là số mũ, nên Prover cam kết
   cột m như một cột witness bình thường. Đổi lại, (3) chỉ tương đương với
   "f ⊂ t" khi đặc số p lớn hơn số lần một giá trị có thể lặp (|f| < p):
   nếu một giá trị v ∉ t lặp đúng p lần, cực điểm -v có hệ số p ≡ 0 và biến
   mất khỏi vế trái (characteristic_counterexample).

2. Bản đơn biến trên miền H = <ω>, |H| = n (cài trong mini-PLONK, Ch.6).
   Cột f (ở đây là cột dây a), selector q_lk (hàng nào phải tra bảng), cột
   bảng t (tiền xử lý, đệm bằng cách lặp phần tử cuối tới n), cột bội m.
   Sau khi có thách thức λ, Prover dựng tổng tích lũy ("running sum") φ:

        φ(ω^0) = 0,  φ(ω^{i+1}) = φ(ω^i) + q_i/(λ + f_i) - m_i/(λ + t_i)

   và chứng minh hai ràng buộc trên H (nhân mẫu để thành đa thức):

        C_1(X) = (φ(ωX) - φ(X))(λ + f(X))(λ + t(X))
                 - q(X)(λ + t(X)) + m(X)(λ + f(X))          = 0 trên H
        C_2(X) = L_0(X) φ(X)                                = 0 trên H

   Vì ràng buộc chuyển tiếp áp cả ở hàng n-1 (ω^n = ω^0), tổng mọi bước bằng
   0, tức (3) đúng tại X = λ. Hai ràng buộc cộng vào tử số của t(X) với hệ số
   α^3, α^4 (lookup thứ k dùng α^{3+2k}, α^{4+2k}).

   Gắn vào mini-PLONK theo api của C06.0 (lớp con LogupProver):
     vòng 1  : thêm [m_k]_1 (m làm mù 2 hệ số)
     vòng 2b : sau [z]_1 rút λ, thêm [φ_k]_1 (φ làm mù 3 hệ số, mở ở ζ và ζω)
     vòng 3  : tử số + α^{3+2k} C_1 + α^{4+2k} C_2
     vòng 4  : thêm m̄, φ̄, φ̄ω, t̄, q̄ (t, q là cột tiền xử lý)
     vòng 5  : W_ζ gộp thêm (m, φ, t, q) với v^6.., W_ζω gộp z và φ_k (v^{k+1})
   Verifier (verify) chép bước 1-12 của plonk.verify, cộng hằng số lookup vào
   r_0 và các cam kết mới vào MSM.
   Proof: 9 + 2K phần tử G1, 6 + 5K phần tử F với K lookup.

3. Kiểm phạm vi (api cho Ch.12): range_check (đồng nhất thức (3) tại λ từ
   transcript), mach_phan_tich_lookup(bounds) (mạch a·b = N với a, b tra
   bảng {0..B}), accepted, vet_can.
"""

import random
from collections import namedtuple

from .field import GF
from .poly import Poly
from .transcript import Transcript
from .curve import R_TOY
from .plonkish import PlonkishCircuit
from . import plonk
from . import opcount as oc
from .multiset import multiplicities, lookup_product_check, product_at

R = R_TOY
F = GF(R)
LABEL = b"snarklib/plonk-logup"


# ---------------------------------------------------------------------------
# 1. Đồng nhất thức LogUp ở mức trường
# ---------------------------------------------------------------------------

def field_multiplicities(f, t, sel=None):
    """m_j ∈ F (bội lấy mod p) cho các f_i có sel_i ≠ 0; trả về (m, ngoai_bang)."""
    fs = [v for i, v in enumerate(f) if sel is None or int(sel[i]) != 0]
    m, outside = multiplicities(fs, t)
    return [F(k) for k in m], outside


def lhs_at(f, X, sel=None):
    """Σ_i sel_i/(X + f_i) (sel_i = 1 nếu sel là None)."""
    X, acc = F(X), F(0)
    for i, v in enumerate(f):
        s = F(1) if sel is None else F(sel[i])
        if s != 0:
            acc = acc + s / (X + F(v))
    return acc


def rhs_at(t, m, X):
    """Σ_j m_j/(X + t_j)."""
    X = F(X)
    return sum((F(k) / (X + F(v)) for v, k in zip(t, m) if F(k) != 0), F(0))


def identity_at(f, t, m, X, sel=None):
    """(3) tại X: True/False; ZeroDivisionError nếu X trùng -f_i hoặc -t_j."""
    return lhs_at(f, X, sel) == rhs_at(t, m, X)


def pole_coefficients(f, t, m, sel=None):
    """Hệ số của 1/(X + v) trong vế trái trừ vế phải, theo từng giá trị v.

    (3) đúng như đồng nhất thức phân thức ⟺ mọi hệ số bằng 0 trong F_p,
    tức count_f(v) ≡ m(v) (mod p). Trả về dict v -> hệ số (chỉ hệ số ≠ 0).
    """
    c = {}
    for i, v in enumerate(f):
        if sel is None or int(sel[i]) != 0:
            c[int(F(v))] = c.get(int(F(v)), F(0)) + 1
    for v, k in zip(t, m):
        c[int(F(v))] = c.get(int(F(v)), F(0)) - F(k)
    return {v: k for v, k in c.items() if k != 0}


def running_sum(f, t, m, lam, sel=None):
    """φ_0..φ_n với φ_{i+1} = φ_i + sel_i/(λ + f_i) - m_i/(λ + t_i).

    f, t, m, sel là các cột cùng độ dài n (một giá trị mỗi hàng của H).
    φ_n = 0 ⟺ (3) đúng tại λ. ZeroDivisionError nếu một mẫu cần dùng bằng 0.
    """
    lam, phi = F(lam), [F(0)]
    for i in range(len(t)):
        s = F(1) if sel is None else F(sel[i])
        dt = lam + F(t[i])
        if dt == 0:
            raise ZeroDivisionError("λ + t_i = 0")
        step = -F(m[i]) / dt
        if s != 0:
            df = lam + F(f[i])
            if df == 0:
                raise ZeroDivisionError("λ + f_i = 0")
            step = step + s / df
        phi.append(phi[-1] + step)
    return phi


def pad_table(t, n):
    """Đệm bảng tới n hàng bằng cách lặp phần tử cuối (bội dồn vào lần đầu)."""
    t = [F(v) for v in t]
    if len(t) > n:
        raise ValueError("bảng dài hơn miền H")
    return t + [t[-1]] * (n - len(t))


def range_table(bound):
    """Bảng {0, 1, ..., bound}."""
    return [F(v) for v in range(bound + 1)]


def range_check(values, bound, label=b"ch07/range", tries=8):
    """Kiểm 'mọi v ∈ values thuộc {0..bound}' bằng đồng nhất thức LogUp (3).

    Prover trung thực tính m từ các giá trị nằm trong bảng; giá trị ngoài bảng
    không có bội tương ứng nên (3) sai tại λ trừ xác suất ≈ (|f|+|t|)/p.
    λ rút từ transcript băm (values, bảng, m); nếu λ chạm cực điểm thì rút lại.
    Trả về (ok, info) với info = {m, outside, lam}.
    """
    t = range_table(bound)
    f = [F(v) for v in values]
    m, outside = field_multiplicities(f, t)
    tr = Transcript(label, R)
    tr.absorb(b"f", f)
    tr.absorb(b"t", t)
    tr.absorb(b"m", m)
    for k in range(tries):
        lam = tr.challenge(b"lambda")
        try:
            ok = identity_at(f, t, m, lam)
        except ZeroDivisionError:
            continue
        return ok, {"m": m, "outside": outside, "lam": lam}
    raise RuntimeError("λ chạm cực điểm quá nhiều lần")


def characteristic_counterexample(table=None, honest=(17, 23), bad=40):
    """Ví dụ (3) 'đúng' dù có giá trị ngoài bảng: bad lặp đúng p = 449 lần.

    f = honest ‖ [bad]·p, bảng {0..31}, m tính từ phần trung thực. Hệ số cực
    điểm của bad là p ≡ 0 nên LogUp chấp nhận tại MỌI λ; kiểm tích (2) (số mũ
    nguyên, UFD) vẫn bác. Trên F_449 cần |f| ≥ p = 449 lookup -- lớn hơn mọi
    miền NTT có được (2-adicity 6, |H| ≤ 64), nên không dựng được trong PIOP.
    """
    table = range_table(31) if table is None else [F(v) for v in table]
    f = [F(v) for v in honest] + [F(bad)] * R
    m, outside = field_multiplicities(f, table)
    mint, _ = multiplicities(f, table)
    rng = random.Random(1530)
    lams = [F(rng.randrange(R)) for _ in range(50)]
    lams = [l for l in lams if l + F(bad) != 0 and all(l + v != 0 for v in table)]
    logup_ok = all(identity_at(f, table, m, l) for l in lams)
    prod_ok = any(lookup_product_check(F, f, table, mint, l) for l in lams)
    return {"f_len": len(f), "outside_count": len(outside), "bad": bad,
            "poles": pole_coefficients(f, table, m), "logup_ok": logup_ok,
            "product_ok_any": prod_ok, "lams": len(lams)}


def min_cheat_size(p=R):
    """Số lookup tối thiểu để giá trị ngoài bảng lọt (3): bội phải ≡ 0 mod p, ≥ p."""
    return p


# ---------------------------------------------------------------------------
# 2. Mạch Plonkish có lookup + mini-PLONK mở rộng
# ---------------------------------------------------------------------------

class LookupCircuit:
    """Mạch Plonkish + danh sách lookup (tên, các hàng tra bảng ở cột a, bảng)."""

    def __init__(self, circ, lookups):
        self.circ = circ
        self.lookups = [(nm, list(rows), [F(v) for v in tab]) for nm, rows, tab in lookups]
        self.n_public = circ.n_public

    def domain(self):
        rows = max([self.circ.n] + [len(t) for _, _, t in self.lookups])
        return plonk.domain_size(rows)


def mach_phan_tich_lookup(bounds=(31, 31)):
    """Mạch a·b = N với kiểm phạm vi bằng lookup thay cho phân tích bit.

    Hàng 0: cổng công khai a = N (q_L = 1, PI = -N)
    Hàng 1: cổng nhân a·b = c (q_M = 1, q_O = -1), sao chép c_1 = a_0
    Hàng 2: không cổng, a_2 = b (sao chép a_2 = b_1)
    Lookup: bounds = (B, B) -> một bảng {0..B} tra ở hàng 1 và 2;
            bounds = (Ba, Bb) khác nhau -> hai lookup: hàng 1 vào {0..Ba},
            hàng 2 vào {0..Bb} (ví dụ (19, 23) của Ch.12).
    """
    rows = [{"qL": 1}, {"qM": 1, "qO": -1}, {}]
    copies = [(("c", 1), ("a", 0)), (("b", 1), ("a", 2))]
    circ = PlonkishCircuit(F, rows, copies, n_public=1)
    Ba, Bb = bounds
    if Ba == Bb:
        lk = [("ab", [1, 2], range_table(Ba))]
    else:
        lk = [("a", [1], range_table(Ba)), ("b", [2], range_table(Bb))]
    return LookupCircuit(circ, lk)


def assignment(N, a, b):
    """Ba cột (a, b, c) và x cho mach_phan_tich_lookup."""
    N, a, b = F(N), F(a), F(b)
    return [N, a, b], [F(0), b, F(0)], [F(0), N, F(0)], [N]


def setup(lc, tau, n=None):
    """(pk, vk) của mini-PLONK + cột t_k, q_k tiền xử lý cho mỗi lookup."""
    n = n or lc.domain()
    pk, vk = plonk.setup(lc.circ, tau, n=n)
    lks, vlks = [], []
    with oc.scope():                                   # tiền xử lý không đếm
        for nm, rows, tab in lc.lookups:
            t = pad_table(tab, n)
            q = [F(1) if i in rows else F(0) for i in range(n)]
            tp = Poly(F, plonk.intt(F, t))
            qp = Poly(F, plonk.intt(F, q))
            cm_t, cm_q = plonk.kzg.commit(pk.srs, tp), plonk.kzg.commit(pk.srs, qp)
            lks.append(dict(name=nm, rows=rows, table_raw=tab, table=t, sel=q,
                            t_poly=tp, q_poly=qp,
                            coset_t=[plonk.eval_on_coset(tp, s, n) for s in pk.shifts],
                            coset_q=[plonk.eval_on_coset(qp, s, n) for s in pk.shifts]))
            vlks.append(dict(name=nm, cm_t=cm_t, cm_q=cm_q))
    pk.lk, vk.lk, pk.lc = lks, vlks, lc               # pk.vk là cùng đối tượng vk
    return pk, vk


def _absorb_pre(tr, vk):
    tr.absorb(b"lookup-preprocessed", [c for L in vk.lk for c in (L["cm_t"], L["cm_q"])])


LkProof = namedtuple("LkProof", ["m", "phi", "m_bar", "phi_bar", "phi_omega_bar",
                                 "t_bar", "q_bar"])
LogupProof = namedtuple("LogupProof", ["base", "lk"])


def proof_size(proof):
    """(số phần tử G1, số phần tử F)."""
    K = len(proof.lk)
    return 9 + 2 * K, 6 + 5 * K


class LogupProver(plonk.Prover):
    """Prover mini-PLONK có thêm K lookup trên cột a (xem docstring đầu module)."""

    def __init__(self, pk, a, b, c, x, rng=None, label=LABEL, hash_public=True,
                 check=True):
        super().__init__(pk, a, b, c, x, rng, label, hash_public, check)
        _absorb_pre(self.tr, pk.vk)
        self.m_vals, self.outside = [], []
        for L in pk.lk:
            m, out = field_multiplicities(self.w[0], L["table"], L["sel"])
            if out and check:
                raise ValueError("giá trị ngoài bảng %s: %s" % (L["name"], out))
            self.m_vals.append(m)
            self.outside.append(out)

    # vòng 1: thêm [m_k]_1 --------------------------------------------------
    def round1(self):
        super().round1()
        rb = lambda: F(self.rng.randrange(R))
        self.m_p = [self._blinded(m, [rb(), rb()]) for m in self.m_vals]
        self.cm_m = [self._commit(p) for p in self.m_p]
        self.tr.absorb(b"round1-lookup", self.cm_m)
        for k, cm in enumerate(self.cm_m):
            self._log(1, "[m_%d]_1" % k, cm)

    # vòng 2b: λ và tổng tích lũy φ_k ----------------------------------------
    def round2(self):
        super().round2()
        n = self.n
        self.lam = self.tr.challenge(b"lambda")
        self._log(2, "lambda", self.lam)
        self.phi_vals, self.phi_p = [], []
        rb = lambda: F(self.rng.randrange(R))
        for L, m in zip(self.pk.lk, self.m_vals):
            phi = running_sum(self.w[0], L["table"], m, self.lam, L["sel"])
            if self.check and phi[n] != 0:
                raise ValueError("tổng LogUp khác 0")
            self.phi_vals.append(phi)
            self.phi_p.append(self._blinded(phi[:n], [rb(), rb(), rb()]))
        self.cm_phi = [self._commit(p) for p in self.phi_p]
        self.tr.absorb(b"round2-lookup", self.cm_phi)
        for k, cm in enumerate(self.cm_phi):
            self._log(2, "[phi_%d]_1" % k, cm)

    # vòng 3: thêm hai ràng buộc mỗi lookup ---------------------------------
    def extra_coset_evals(self, j, s):
        out, n = {}, self.n
        for k, L in enumerate(self.pk.lk):
            ph = plonk.eval_on_coset(self.phi_p[k], s, n)
            out["phi%d" % k], out["phiw%d" % k] = ph, ph[1:] + ph[:1]
            out["m%d" % k] = plonk.eval_on_coset(self.m_p[k], s, n)
            out["t%d" % k], out["q%d" % k] = L["coset_t"][j], L["coset_q"][j]
        return out

    def numerator_on_coset(self, j, ev):
        out, lam, al = super().numerator_on_coset(j, ev), self.lam, self.alpha
        for k in range(len(self.pk.lk)):
            a3, a4 = al ** (3 + 2 * k), al ** (4 + 2 * k)
            for i in range(self.n):
                A, T = lam + ev["a"][i], lam + ev["t%d" % k][i]
                ph = ev["phi%d" % k][i]
                C1 = ((ev["phiw%d" % k][i] - ph) * A * T - ev["q%d" % k][i] * T
                      + ev["m%d" % k][i] * A)
                out[i] = out[i] + a3 * C1 + a4 * ev["L0"][i] * ph
        return out

    def quotient_naive(self):
        q, rem = super().quotient_naive()
        num, lam, al = q * self.ZH + rem, self.lam, self.alpha
        for k, L in enumerate(self.pk.lk):
            ph, A, T = self.phi_p[k], self.a_p + lam, L["t_poly"] + lam
            C1 = (ph.scale(self.pk.omega) - ph) * A * T - L["q_poly"] * T + self.m_p[k] * A
            num = num + C1 * al ** (3 + 2 * k) + (self.pk.polys["L0"] * ph) * al ** (4 + 2 * k)
        return num.divmod(self.ZH)

    # vòng 4: thêm 5 giá trị mỗi lookup -------------------------------------
    def round4(self):
        super().round4()
        z = self.zeta
        zw = z * self.pk.omega
        self.lk_evals = [(self.m_p[k](z), self.phi_p[k](z), self.phi_p[k](zw),
                          L["t_poly"](z), L["q_poly"](z)) for k, L in enumerate(self.pk.lk)]
        self.tr.absorb(b"round4-lookup", [v for e in self.lk_evals for v in e])
        for k, e in enumerate(self.lk_evals):
            for nm, v in zip(LkProof._fields[2:], e):
                self._log(4, "%s_%d" % (nm, k), v)

    def lookup_constant(self):
        """Σ_k α^{3+2k} C_1(ζ) + α^{4+2k} L_0(ζ) φ̄_k (hằng số cộng vào r_0)."""
        return lookup_constant(self.lk_evals, self.evals[0], self.lam, self.alpha,
                               plonk.lagrange0_at(self.n, self.zeta))

    # vòng 5 ----------------------------------------------------------------
    def linearization(self):
        return super().linearization() + self.lookup_constant()

    def zeta_openings(self):
        out = super().zeta_openings()
        for k, L in enumerate(self.pk.lk):
            mb, pb, _, tb, qb = self.lk_evals[k]
            out += [(self.m_p[k], mb), (self.phi_p[k], pb), (L["t_poly"], tb),
                    (L["q_poly"], qb)]
        return out

    def round5(self):
        z, om = self.zeta, self.pk.omega
        self.v = self.tr.challenge(b"v")
        self._log(5, "v", self.v)
        acc, vp = Poly.zero(F), F(1)
        for f, y in self.zeta_openings():
            acc = acc + (f - y) * vp
            vp = vp * self.v
        Wz, rem = acc.divmod(Poly(F, [-z, 1]))
        accw, vp = self.z_p - self.evals[5], self.v
        for k in range(len(self.pk.lk)):                # φ_k mở ở ζω, hệ số v^{k+1}
            accw = accw + (self.phi_p[k] - self.lk_evals[k][2]) * vp
            vp = vp * self.v
        Wzw, rem2 = accw.divmod(Poly(F, [-z * om, 1]))
        if self.check:
            assert rem.is_zero() and rem2.is_zero(), "r(ζ) ≠ 0"
        self.cm_W = [self._commit(Wz), self._commit(Wzw)]
        self.tr.absorb(b"round5", self.cm_W)
        self._log(5, "[W_zeta]_1", self.cm_W[0])
        self._log(5, "[W_zeta_omega]_1", self.cm_W[1])
        self.u = self.tr.challenge(b"u")
        self._log(5, "u", self.u)

    def proof(self):
        lk = [LkProof(self.cm_m[k], self.cm_phi[k], *self.lk_evals[k])
              for k in range(len(self.pk.lk))]
        return LogupProof(super().proof(), lk)


def lookup_constant(lk_evals, a_bar, lam, alpha, l0):
    out = F(0)
    for k, (mb, pb, pwb, tb, qb) in enumerate(lk_evals):
        A, T = lam + a_bar, lam + tb
        C1 = (pwb - pb) * A * T - qb * T + mb * A
        out = out + alpha ** (3 + 2 * k) * C1 + alpha ** (4 + 2 * k) * l0 * pb
    return out


def prove(pk, a, b, c, x, rng=None, label=LABEL, hash_public=True, check=True):
    """Trả về (LogupProof, log). check=False: Prover gian (bỏ phần dư của t)."""
    p = LogupProver(pk, a, b, c, x, rng, label, hash_public, check)
    pf = p.run()
    return pf, p.log


def challenges(vk, x, proof, label=LABEL, hash_public=True):
    """β, γ, λ, α, ζ, v, u -- cùng thứ tự transcript với LogupProver."""
    P = proof.base
    tr = plonk.new_transcript(vk, x, label, hash_public)
    _absorb_pre(tr, vk)
    tr.absorb(b"round1", [P.a, P.b, P.c])
    tr.absorb(b"round1-lookup", [L.m for L in proof.lk])
    be, ga = tr.challenge(b"beta"), tr.challenge(b"gamma")
    tr.absorb(b"round2", P.z)
    lam = tr.challenge(b"lambda")
    tr.absorb(b"round2-lookup", [L.phi for L in proof.lk])
    al = tr.challenge(b"alpha")
    tr.absorb(b"round3", [P.t_lo, P.t_mid, P.t_hi])
    ze = tr.challenge(b"zeta")
    tr.absorb(b"round4", [getattr(P, k) for k in plonk.FR_FIELDS])
    tr.absorb(b"round4-lookup", [F(int(getattr(L, k))) for L in proof.lk
                                  for k in LkProof._fields[2:]])
    v = tr.challenge(b"v")
    tr.absorb(b"round5", [P.W_zeta, P.W_zeta_omega])
    u = tr.challenge(b"u")
    return be, ga, lam, al, ze, v, u


def verify(vk, x, proof, label=LABEL, hash_public=True, trace=None):
    """Verifier mini-PLONK + LogUp: bước 1-12 của plonk.verify với phần lookup."""
    P, n, om, k1, k2, cm = proof.base, vk.n, vk.omega, vk.k1, vk.k2, vk.comms
    K = len(vk.lk)
    # 1-3. Kiểm 9 + 2K phần tử G1, 6 + 5K phần tử F, đầu vào công khai
    if len(proof.lk) != K:
        return False
    g1s = [getattr(P, k) for k in plonk.G1_FIELDS] + [Q for L in proof.lk for Q in (L.m, L.phi)]
    if not all(plonk._in_group(Q) for Q in g1s):
        return False
    a, b, c, s1, s2, zw = (F(int(getattr(P, k))) for k in plonk.FR_FIELDS)
    lke = [tuple(F(int(getattr(L, k))) for k in LkProof._fields[2:]) for L in proof.lk]
    if len(x) != vk.ell:
        return False
    x = [F(v) for v in x]
    # 4. Thách thức
    be, ga, lam, al, ze, v, u = challenges(vk, x, proof, label, hash_public)
    # 5-7.
    zh = ze ** n - 1
    l0 = plonk.lagrange0_at(n, ze)
    pi = plonk.pi_at(n, om, x, ze)
    # 8. r_0 = r_0 của PLONK + hằng số lookup (C_1, C_2 tính trọn tại ζ)
    r0 = pi - l0 * al * al - al * (a + be * s1 + ga) * (b + be * s2 + ga) * (c + ga) * zw
    r0 = r0 + lookup_constant(lke, a, lam, al, l0)
    # 9. [D]_1 như PLONK
    sc_z = (a + be * ze + ga) * (b + be * k1 * ze + ga) * (c + be * k2 * ze + ga) * al \
        + l0 * al * al + u
    sc_s3 = -(a + be * s1 + ga) * (b + be * s2 + ga) * al * be * zw
    D = [(cm["qM"], a * b), (cm["qL"], a), (cm["qR"], b), (cm["qO"], c), (cm["qC"], F(1)),
         (P.z, sc_z), (cm["s3"], sc_s3), (P.t_lo, -zh),
         (P.t_mid, -zh * ze ** n), (P.t_hi, -zh * ze ** (2 * n))]
    # 10. [F]_1 = [D]_1 + v[a] + ... + v^5[s2] + Σ_k (v^{6+4k}[m_k] + v^{7+4k}[φ_k]
    #             + v^{8+4k}[t_k] + v^{9+4k}[q_k]) + Σ_k u·v^{k+1}[φ_k] (mở ở ζω)
    Fl = D + [(P.a, v), (P.b, v ** 2), (P.c, v ** 3), (cm["s1"], v ** 4), (cm["s2"], v ** 5)]
    # 11. [E]_1
    e = -r0 + v * a + v ** 2 * b + v ** 3 * c + v ** 4 * s1 + v ** 5 * s2 + u * zw
    for k, L in enumerate(proof.lk):
        mb, pb, pwb, tb, qb = lke[k]
        vb = v ** (6 + 4 * k)
        Fl += [(L.m, vb), (L.phi, vb * v + u * v ** (k + 1)),
               (vk.lk[k]["cm_t"], vb * v ** 2), (vk.lk[k]["cm_q"], vb * v ** 3)]
        e = e + vb * (mb + v * pb + v ** 2 * tb + v ** 3 * qb) + u * v ** (k + 1) * pwb
    # 12. e([W_ζ] + u[W_ζω], [x]_2) = e(ζ[W_ζ] + uζω[W_ζω] + [F] - [E], [1]_2)
    pts = [Q for Q, _ in Fl] + [vk.G1, P.W_zeta, P.W_zeta_omega]
    scs = [s for _, s in Fl] + [-e, ze, u * ze * om]
    rhs = oc.msm(pts, scs)
    lhs = oc.msm([P.W_zeta, P.W_zeta_omega], [F(1), u])
    ok = oc.pairing_check([(lhs, vk.g2[1]), (-rhs, vk.g2[0])])
    if trace is not None:
        trace.update(dict(beta=be, gamma=ga, lam=lam, alpha=al, zeta=ze, v=v, u=u,
                          ZH=zh, L0=l0, PI=pi, r0=r0, E=e, ok=ok, msm_size=len(pts)))
    return ok


# ---------------------------------------------------------------------------
# 3. API kiểm phạm vi cho Ch.12
# ---------------------------------------------------------------------------

def accepted(N, a, b, bounds=(31, 31)):
    """Witness (a, b) có thỏa mạch mach_phan_tich_lookup(bounds) không (không cần proof):
    cổng a·b = N trong F_449 và a ∈ {0..Ba}, b ∈ {0..Bb}."""
    a, b = F(a), F(b)
    return a * b == F(N) and int(a) <= bounds[0] and int(b) <= bounds[1]


def vet_can(N, bounds=(31, 31)):
    """Mọi (a, b) ∈ F_449^2 thỏa mạch -- vét a trong bảng, b = N/a."""
    out = []
    for a in range(1, bounds[0] + 1):
        b = int(F(N) / F(a))
        if b <= bounds[1]:
            out.append((a, b))
    return out


def prove_factor(N, a, b, bounds=(31, 31), tau=123, rng=None, label=LABEL, check=True):
    """Tiện ích một dòng: setup + prove + verify cho mạch a·b = N có lookup.

    Trả về (ok, proof, pk, vk). check=True và witness sai -> ValueError.
    """
    lc = mach_phan_tich_lookup(bounds)
    pk, vk = setup(lc, tau)
    A, B, C, x = assignment(N, a, b)
    pf, _ = prove(pk, A, B, C, x, rng or random.Random(0), label=label, check=check)
    return verify(vk, x, pf, label=label), pf, pk, vk


def circuit_stats(lc):
    """Số hàng cổng, số lookup, cỡ bảng, n, số cam kết tiền xử lý, cỡ proof."""
    K = len(lc.lookups)
    return {"gate_rows": lc.circ.n, "lookups": sum(len(r) for _, r, _ in lc.lookups),
            "tables": [len(t) for _, _, t in lc.lookups], "n": lc.domain(),
            "preprocessed_comms": 8 + 2 * K, "proof_G1": 9 + 2 * K, "proof_F": 6 + 5 * K,
            "K": K}
