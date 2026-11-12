"""Chương 9 — Mã phụ (không in): tích lũy phép kiểm IPA kiểu Halo (Mục 9.2b, eq:ch9-halo-acc).

Dùng PCS từ IPA của snarklib/ipa.py (b công khai, không blinding) trên nhóm cấp 449.
(1) Verifier IPA phải tính G_f = <s, G> (MSM cỡ n) và b_f = <s, (1, z, ..., z^{n-1})>.
    s là hệ số của g(X) = Π_j (u_j^{-1} + u_j X^{2^{k-j}}), nên b_f = g(z) tính trong
    O(log n), còn G_f = Commit(g) -- phần duy nhất tốn O(n).
(2) Bộ tích lũy acc = (G_f, u_1..u_k). Gộp hai bộ tích lũy: rút alpha, x; khẳng định
    Commit(g_1 + alpha g_2) = G_1 + alpha G_2 mở tại x ra g_1(x) + alpha g_2(x); kiểm một
    lần mở IPA thay cho hai MSM cỡ n. Kiểm âm: G_2 sai (cộng thêm một điểm sinh).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.transcript import Transcript
from snarklib.pedersen import msm, smul, padd
from snarklib import ipa

F = ipa.F
n, k = 8, 3
pp = ipa.setup(n, b"ch9-halo")


def challenges(P, z, v, proof, label):
    """Tái tạo các thách thức u_1..u_k như Verifier (Fiat--Shamir)."""
    t = Transcript(label)
    ipa._bind(t, pp, P, z, v)
    us = []
    for L, R in zip(proof["L"], proof["R"]):
        t.absorb(b"L", L)
        t.absorb(b"R", R)
        us.append(ipa.challenge_nonzero(t, b"u"))
    return us


def poly_g(us):
    """Hệ số của Π_j (u_j^{-1} + u_j X^{2^{k-j}}) (j = 1..k), bậc n - 1."""
    c = [F(1)]
    for j, u in enumerate(us):
        e = 1 << (len(us) - 1 - j)
        d = [F(0)] * (len(c) + e)
        for i, ci in enumerate(c):
            d[i] += ci * u.inverse()
            d[i + e] += ci * u
        c = d
    return c


def mot_proof(a, z, label):
    P = ipa.commit(pp, a)
    v, pr = ipa.open(pp, a, z, Transcript(label))
    assert ipa.verify(pp, P, z, v, pr, Transcript(label))
    us = challenges(P, z, v, pr, label)
    s = ipa.s_vector(us)
    assert poly_g(us) == s                                   # s = hệ số của g(X)
    assert ipa.inner(s, ipa.powers(z, n)) == ipa.b_final_powers(us, z)
    return us, s, msm(s, pp.G)


a1 = [17, 23, 28, 30, 2, 3, 5, 7]
a2 = [391, 1, 4, 9, 16, 25, 36, 49]
us1, s1, G1 = mot_proof(a1, 5, b"ch9-acc-1")
us2, s2, G2 = mot_proof(a2, 11, b"ch9-acc-2")
print("n = %d, k = %d vòng IPA" % (n, k))
print("proof 1: u =", [int(u) for u in us1], " s = hệ số g_1 =", [int(x) for x in s1])
print("proof 2: u =", [int(u) for u in us2], " s = hệ số g_2 =", [int(x) for x in s2])
print("G_1 = <s_1, G> =", G1, " G_2 =", G2)
print("b_f = g(z): proof 1 g_1(5) = %d, proof 2 g_2(11) = %d (tích %d thừa số thay vì %d số hạng)"
      % (ipa.b_final_powers(us1, 5), ipa.b_final_powers(us2, 11), k, n))


def tich_luy(Ga, usa, Gb, usb, label):
    """Gộp hai bộ tích lũy thành một khẳng định mở; trả (alpha, x, C, v, kết quả kiểm)."""
    t = Transcript(label)
    for P in (Ga, Gb):
        t.absorb(b"G", P)
    alpha, x = t.challenge(b"alpha"), t.challenge(b"x")
    C = padd(Ga, smul(Gb, alpha))                            # Verifier: 1 nhân vô hướng
    v = ipa.b_final_powers(usa, x) + alpha * ipa.b_final_powers(usb, x)   # O(log n)
    coef = [p + alpha * q for p, q in zip(poly_g(usa), poly_g(usb))]     # helper / Prover
    v2, pr = ipa.open(pp, coef, x, Transcript(label + b"/open"))
    cnt = ipa.VCount()
    ok = (v2 == v) and ipa.verify(pp, C, x, v, pr, Transcript(label + b"/open"), cnt=cnt)
    return alpha, x, C, v, ok, cnt


alpha, x, C, v, ok, cnt = tich_luy(G1, us1, G2, us2, b"ch9-acc-gop")
print("gộp: alpha = %d, x = %d, C = G_1 + alpha G_2 = %s, v = g_1(x) + alpha g_2(x) = %d"
      % (alpha, x, C, v))
print("mở IPA cho C tại x: hợp lệ =", ok, "| Verifier:", cnt.as_dict())
print("so sánh: kiểm riêng hai bộ tích lũy = 2 MSM cỡ %d (%d nhân vô hướng); gộp = 1 lần mở" % (n, 2 * n))
assert ok

# Kiểm âm: bộ tích lũy thứ hai sai (G_2 lệch một điểm sinh), thử 20 nhãn transcript
G2_sai = padd(G2, pp.G[0])
lot = sum(tich_luy(G1, us1, G2_sai, us2, b"ch9-acc-sai-%d" % i)[4] for i in range(20))
print("G_2 sai: lọt %d/20 lần gộp (trên trường 449 xác suất lọt mỗi lần cỡ 1/449)" % lot)
