"""PCS kiểu Hyrax cho đa thức đa tuyến tính: cam kết ma trận √n × √n.

Bảng giá trị T (n = 2^l phần tử, thứ tự big-endian của snarklib.mle) được xếp
thành ma trận M cỡ m × k (m = 2^{l_r}, k = 2^{l_c}, l_r = l // 2, l_c = l - l_r):
    M[i][j] = T[i·k + j]  — chỉ số hàng i ứng với các biến CAO x_1..x_{l_r}.

Với điểm r = (r_hàng, r_cột):
    f~(r) = Σ_{i,j} M[i][j]·eq(i, r_hàng)·eq(j, r_cột) = L^T · M · R,
    L = eq_table(r_hàng) (độ dài m),   R = eq_table(r_cột) (độ dài k).

Cam kết: C_i = <M[i], G> cho từng hàng (m điểm, không blinding ở bản này).
Mở:  Verifier tự tính C_w = Σ_i L_i·C_i = <w, G> với w = L^T·M (tính đồng cấu),
     Prover chứng minh <w, R> = v bằng IPA với b = R công khai (snarklib.ipa),
     hoặc (chế độ "truc_tiep") gửi thẳng w (k phần tử trường).
Kích thước: commitment m điểm; proof IPA 2·log2 k điểm + 1 vô hướng.

Hyrax gốc (Wahby et al. 2018) có blinding và "proof-of-dot-product" kiểu
Bulletproofs để đạt zero-knowledge; bản ở đây giữ phần cấu trúc √n, lược bỏ
blinding (ghi rõ trong sách).
"""

from .field import GF
from .curve import R_TOY
from .mle import eq_table, num_vars
from .pedersen import PedersenParams, msm, smul, padd
from . import ipa

F = GF(R_TOY)


def shape(n):
    """(l_r, l_c, m, k) cho bảng n phần tử."""
    l = num_vars([0] * n)
    l_r = l // 2
    l_c = l - l_r
    return l_r, l_c, 1 << l_r, 1 << l_c


def setup(n, label=b"snarklib/hyrax"):
    _, _, _, k = shape(n)
    return PedersenParams(k, label)


def to_matrix(table):
    _, _, m, k = shape(len(table))
    return [[F(table[i * k + j]) for j in range(k)] for i in range(m)]


def commit(pp, table):
    """Trả về (danh sách m cam kết hàng, trạng thái của Prover)."""
    M = to_matrix(table)
    comms = [msm(row, pp.G[:len(row)]) for row in M]
    return comms, {"M": M, "n": len(table)}


def _split(n, r):
    l_r, l_c, _, _ = shape(n)
    if len(r) != l_r + l_c:
        raise ValueError("số biến không khớp")
    return eq_table(F, r[:l_r]), eq_table(F, r[l_r:])


def combine_rows(comms, L, cnt=None):
    """C_w = Σ_i L_i·C_i (Verifier tự tính nhờ tính đồng cấu)."""
    return msm(L, comms, cnt)


def open(pp, state, r, transcript, mode="ipa"):
    """Mở tại r. Trả về (v, proof)."""
    M, n = state["M"], state["n"]
    L, R = _split(n, r)
    k = len(R)
    w = [sum((L[i] * M[i][j] for i in range(len(M))), F(0)) for j in range(k)]
    v = ipa.inner(w, R)
    transcript.absorb(b"r", [F(x) for x in r])
    transcript.absorb(b"v", v)
    if mode == "truc_tiep":
        return v, {"w": w}
    xu = ipa.challenge_nonzero(transcript, b"x_u")
    U1 = smul(pp.U, xu)
    return v, {"ipa": ipa.prove_public_b(pp.G[:k], U1, w, R, transcript)}


def verify(pp, comms, r, v, proof, transcript, n, cnt=None):
    L, R = _split(n, r)
    k = len(R)
    if len(comms) != len(L):
        return False
    C_w = combine_rows(comms, L, cnt)
    transcript.absorb(b"r", [F(x) for x in r])
    transcript.absorb(b"v", F(v))
    if "w" in proof:                                   # chế độ trực tiếp
        w = proof["w"]
        return len(w) == k and msm(w, pp.G[:k], cnt) == C_w and ipa.inner(w, R) == F(v)
    xu = ipa.challenge_nonzero(transcript, b"x_u")
    U1 = smul(pp.U, xu, cnt)
    P1 = padd(C_w, smul(U1, v, cnt), cnt)               # P' = C_w + v·U'
    return ipa.verify_public_b(pp.G[:k], U1, P1, R, proof["ipa"], transcript, cnt)


def sizes(n, mode="ipa"):
    """Kích thước (số phần tử) của commitment và proof."""
    _, l_c, m, k = shape(n)
    if mode == "truc_tiep":
        return {"commit_group": m, "proof_group": 0, "proof_field": k}
    return {"commit_group": m, "proof_group": 2 * l_c, "proof_field": 1}
