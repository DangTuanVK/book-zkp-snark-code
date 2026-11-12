"""Kiểm đa tập bằng tích (Chương 7, §7.1): nền của mọi lookup kiểu tích.

Đa tập (multiset) là tập cho phép phần tử lặp; hai đa tập bằng nhau khi cùng
phần tử với cùng bội (chen2023hyperplonk.txt dòng 456-457). Với hai dãy
f = (f_1..f_n), g = (g_1..g_n) trên F_p:

    {f_i} = {g_i} (đa tập)  ⟺  Π_i (X + f_i) = Π_i (X + g_i)  trong F_p[X]   (1)

vì F_p[X] là miền phân tích duy nhất: hai vế là tích các nhân tử bậc một và
nghiệm -f_i đếm đúng bội (cùng lập luận UFD ở Claim 3.1 của Plookup,
gabizon2020plookup.txt dòng 164-197). Verifier không kiểm (1) như đa thức mà
thay X bằng thách thức ngẫu nhiên γ (Schwartz–Zippel: lọt với xác suất
≤ n/p) -- đúng PIOP "Multiset Check" của HyperPlonk §3.4: f' = r + f,
g' = r + g rồi kiểm tích (chen2023hyperplonk.txt dòng 1095-1118).

Lookup dạng tích: f ⊂ t với bội m_j (số lần t_j xuất hiện trong f) ⟺

    Π_i (X + f_i) = Π_j (X + t_j)^{m_j}                                        (2)

(eq:ch7-multiset). Bội m_j là SỐ MŨ nguyên: vế phải có bậc Σ m_j = n và không
viết được thành ràng buộc bậc thấp theo từng hàng khi m_j tùy ý; vì vậy Plookup
phải sắp xếp (vector s) và LogUp chuyển sang tổng phân thức (snarklib.logup),
nơi m_j trở thành HỆ SỐ trong trường.

Khác LogUp: (1) và (2) đúng với mọi đặc số p, kể cả khi một giá trị lặp ≥ p
lần (UFD không quan tâm đặc số). Xem logup.characteristic_counterexample.
"""

import random

from .field import GF
from .poly import Poly

F449 = GF(449)


def _F(F, xs):
    return [F(x) for x in xs]


def multiplicities(f, t):
    """Bội m_j = #{i : f_i = t_j} (số nguyên); giá trị lặp trong t chỉ nhận ở
    lần xuất hiện đầu. Trả về (m, ngoai_bang) với ngoai_bang là các f_i ∉ t."""
    pos, m, outside = {}, [0] * len(t), []
    for j, v in enumerate(t):
        pos.setdefault(int(v), j)
    for v in f:
        j = pos.get(int(v))
        if j is None:
            outside.append(int(v))
        else:
            m[j] += 1
    return m, outside


def expand(t, m):
    """Đa tập t^m dưới dạng dãy: t_j lặp m_j lần (m_j nguyên ≥ 0)."""
    out = []
    for v, k in zip(t, m):
        out.extend([v] * int(k))
    return out


def product_at(F, vals, X, exps=None):
    """Π_i (X + v_i)^{e_i} tại một điểm X (e_i = 1 nếu exps là None)."""
    acc, X = F(1), F(X)
    for i, v in enumerate(vals):
        e = 1 if exps is None else int(exps[i])
        acc = acc * (X + F(v)) ** e
    return acc


def product_poly(F, vals, exps=None):
    """Π_i (X + v_i)^{e_i} dưới dạng đa thức (để kiểm đồng nhất thức hình thức)."""
    acc = Poly.const(F, 1)
    for i, v in enumerate(vals):
        e = 1 if exps is None else int(exps[i])
        for _ in range(e):
            acc = acc * Poly(F, [F(v), 1])
    return acc


def multiset_equal_poly(F, f, g):
    """Kiểm (1) như đồng nhất thức đa thức (tốn O(n^2), chỉ để đối chiếu)."""
    return product_poly(F, f) == product_poly(F, g)


def multiset_check(F, f, g, gamma):
    """Kiểm (1) tại X = γ: Π (γ + f_i) == Π (γ + g_i)."""
    return product_at(F, f, gamma) == product_at(F, g, gamma)


def lookup_product_check(F, f, t, m, gamma):
    """Kiểm (2) tại X = γ: Π (γ + f_i) == Π (γ + t_j)^{m_j} (m_j nguyên)."""
    return product_at(F, f, gamma) == product_at(F, t, gamma, m)


def grand_product(F, f, g, gamma):
    """Tích tích lũy Z_0 = 1, Z_{i+1} = Z_i·(γ + f_i)/(γ + g_i), i = 0..n-1.

    Trả về Z_0..Z_n; Z_n = 1 ⟺ (1) đúng tại γ (cùng cơ chế với z(X) của PLONK).
    ZeroDivisionError nếu γ + g_i = 0 (Prover thật rút γ khác).
    """
    Z = [F(1)]
    for a, b in zip(f, g):
        den = F(gamma) + F(b)
        if den == 0:
            raise ZeroDivisionError("γ + g_i = 0")
        Z.append(Z[-1] * (F(gamma) + F(a)) / den)
    return Z


def soundness_bound(n, p):
    """Cận Schwartz–Zippel cho kiểm (1)/(2) tại một γ: deg / p = n / p."""
    return n / p


def false_positive_rate(F, f, g, trials=None, rng=None):
    """Tỉ lệ γ ∈ F_p làm (1) 'đúng' dù f ≠ g. trials=None: vét toàn bộ F_p."""
    p = F.p
    if trials is None:
        gammas = range(p)
    else:
        rng = rng or random.Random(0)
        gammas = [rng.randrange(p) for _ in range(trials)]
    hits = sum(1 for gm in gammas if multiset_check(F, f, g, gm))
    return hits, len(gammas)
