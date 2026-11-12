"""PCS kiểu Ligero cho đa thức đa tuyến tính (mã Reed–Solomon + Merkle).

Bảng giá trị T (n = 2^l, big-endian như snarklib.mle) xếp thành ma trận M cỡ
m × k (m = 2^{l//2} hàng ứng với các biến cao, k = 2^{l - l//2} cột), giống
hyrax.py. Khi đó f~(r) = L^T·M·R với L = eq_table(r_hàng), R = eq_table(r_cột).

Commit:
    mỗi hàng M[i] (k phần tử) xem là hệ số của đa thức bậc < k, đánh giá trên
    miền con nhân H' cỡ ρ^{-1}·k của F_449 (ρ = 1/4 mặc định, dùng snarklib.fft)
    → hàng mã Enc(M[i]) dài ρ^{-1}·k. Cam kết = gốc Merkle trên CÁC CỘT của ma
    trận mã (mỗi lá là một cột m phần tử).
Open(r) — sau khi hấp thụ (gốc, r, v) vào transcript:
    1. kiểm proximity: γ ∈ F^m ngẫu nhiên, Prover gửi u_prox = γ^T·M  (k phần tử);
    2. đánh giá:       Prover gửi u_eval = L^T·M                      (k phần tử);
    3. t chỉ số cột ngẫu nhiên (từ transcript); Prover mở các cột đó + đường Merkle.
Verify:
    <u_eval, R> = v;  với mỗi cột j được mở: đường Merkle đúng,
    Enc(u_prox)[j] = Σ_i γ_i·M'[i][j]  và  Enc(u_eval)[j] = Σ_i L_i·M'[i][j].
Tính tuyến tính của mã cho Enc(γ^T M) = γ^T Enc(M), nên Prover trung thực luôn qua.

Kích thước proof: 2k + t·m phần tử trường, t·log2(ρ^{-1} k) giá trị băm
(không gộp các nút Merkle trùng nhau).

LƯU Ý: trên F_449 thách thức chỉ có 449 giá trị và mã ngắn, nên độ an toàn là
vô nghĩa; mã chỉ minh họa cơ chế. Tham số t thực tế tính theo mức an toàn
(Mục 4.2c). Miền RS lớn nhất trên F_449 là 64 phần tử → k ≤ 64/ρ^{-1}.
"""

from .field import GF
from .curve import R_TOY
from .fft import ntt
from .mle import eq_table, num_vars
from . import merkle

F = GF(R_TOY)
T_DEFAULT = 8
RHO_INV_DEFAULT = 4
HASH_BYTES = 32


def shape(n):
    l = num_vars([0] * n)
    l_r = l // 2
    return l_r, l - l_r, 1 << l_r, 1 << (l - l_r)


def rs_encode(row, rho_inv=RHO_INV_DEFAULT):
    """Mã RS tỉ lệ 1/rho_inv: đánh giá đa thức hệ số `row` trên <ω>, |<ω>| = rho_inv·k."""
    k = len(row)
    return ntt(F, [F(x) for x in row] + [F(0)] * ((rho_inv - 1) * k))


def column_bytes(col):
    w = (F.p.bit_length() + 7) // 8
    return b"".join(int(x).to_bytes(w, "big") for x in col)


def _lin(coefs, col):
    return sum((c * x for c, x in zip(coefs, col)), F(0))


def commit(table, rho_inv=RHO_INV_DEFAULT):
    """Trả về (gốc Merkle, trạng thái của Prover)."""
    n = len(table)
    _, _, m, k = shape(n)
    M = [[F(table[i * k + j]) for j in range(k)] for i in range(m)]
    enc = [rs_encode(row, rho_inv) for row in M]
    ncols = rho_inv * k
    cols = [[enc[i][j] for i in range(m)] for j in range(ncols)]
    root, tree = merkle.commit([column_bytes(c) for c in cols])
    return root, {"M": M, "enc": enc, "cols": cols, "tree": tree,
                  "n": n, "m": m, "k": k, "rho_inv": rho_inv}


def query_indices(transcript, ncols, t):
    """t chỉ số cột KHÁC NHAU rút từ transcript (t bị chặn bởi số cột)."""
    t = min(t, ncols)
    out, ctr = [], 0
    while len(out) < t:
        raw = transcript.challenge_bytes(b"q%d" % ctr, 8)
        j = int.from_bytes(raw, "big") % ncols
        if j not in out:
            out.append(j)
        ctr += 1
    return out


def _prelude(transcript, root, r, v, m):
    transcript.absorb(b"root", root)
    transcript.absorb(b"r", [F(x) for x in r])
    transcript.absorb(b"v", F(v))
    return [transcript.challenge(b"gamma%d" % i) for i in range(m)]


def evaluate(state, r):
    """v = f~(r) = L^T·M·R tính từ trạng thái."""
    l_r = shape(state["n"])[0]
    L, R = eq_table(F, r[:l_r]), eq_table(F, r[l_r:])
    return sum((L[i] * R[j] * state["M"][i][j]
                for i in range(state["m"]) for j in range(state["k"])), F(0))


def open(state, r, transcript, t=T_DEFAULT):
    """Proof mở f~ tại r (giá trị v = evaluate(state, r) cũng được ghi vào proof)."""
    m, k, M = state["m"], state["k"], state["M"]
    l_r = shape(state["n"])[0]
    L = eq_table(F, r[:l_r])
    v = evaluate(state, r)
    gamma = _prelude(transcript, state["tree"].root, r, v, m)
    u_prox = [_lin(gamma, [M[i][j] for i in range(m)]) for j in range(k)]
    u_eval = [_lin(L, [M[i][j] for i in range(m)]) for j in range(k)]
    transcript.absorb(b"u_prox", u_prox)
    transcript.absorb(b"u_eval", u_eval)
    idx = query_indices(transcript, state["rho_inv"] * k, t)
    cols = [(j, state["cols"][j], state["tree"].open(j)) for j in idx]
    return {"v": v, "u_prox": u_prox, "u_eval": u_eval, "cols": cols}


def verify(root, r, value, proof, transcript, t=T_DEFAULT, rho_inv=RHO_INV_DEFAULT, cnt=None):
    """Kiểm proof. cnt (dict, tùy chọn) nhận số phép băm và nhân trường của Verifier."""
    n = 1 << len(r)
    l_r, _, m, k = shape(n)
    L, R = eq_table(F, r[:l_r]), eq_table(F, r[l_r:])
    ncols = rho_inv * k
    u_prox, u_eval = proof["u_prox"], proof["u_eval"]
    if len(u_prox) != k or len(u_eval) != k:
        return False
    gamma = _prelude(transcript, root, r, value, m)
    if sum((a * b for a, b in zip(u_eval, R)), F(0)) != F(value):
        return False
    transcript.absorb(b"u_prox", u_prox)
    transcript.absorb(b"u_eval", u_eval)
    idx = query_indices(transcript, ncols, t)
    e_prox, e_eval = rs_encode(u_prox, rho_inv), rs_encode(u_eval, rho_inv)
    if len(proof["cols"]) != len(idx):
        return False
    hashes = 0
    for j, (jj, col, path) in zip(idx, proof["cols"]):
        if jj != j or len(col) != m:
            return False
        if not merkle.verify(root, j, column_bytes(col), path):
            return False
        hashes += 1 + len(path)
        if _lin(gamma, col) != e_prox[j] or _lin(L, col) != e_eval[j]:
            return False
    if cnt is not None:
        lg = ncols.bit_length() - 1
        cnt["hash"] = hashes
        # 2 NTT cỡ ncols (≈ (ncols/2)·log2(ncols) nhân mỗi cái) + 2m nhân mỗi cột + k cho <u_eval,R>
        cnt["fmul"] = 2 * (ncols // 2) * lg + 2 * m * len(idx) + k
    return True


def sizes(n, t=T_DEFAULT, rho_inv=RHO_INV_DEFAULT):
    """Số phần tử trường và số giá trị băm trong proof; commitment = 1 băm."""
    _, _, m, k = shape(n)
    ncols = rho_inv * k
    t = min(t, ncols)
    return {"m": m, "k": k, "ncols": ncols, "t": t,
            "proof_field": 2 * k + t * m,
            "proof_hash": t * (ncols.bit_length() - 1),
            "commit_hash": 1}
