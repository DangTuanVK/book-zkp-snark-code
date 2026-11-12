"""Mini-Spartan: zkSNARK (không ZK) cho R1CS bằng hai sum-check + PCS đa tuyến tính.

Đối chiếu setty2020spartan.txt:
  - z và R1CS dưới dạng hàm trên siêu lập phương, F_io(x) (:818-826);
  - Q_io(τ) = Σ_x F~_io(x)·eq(τ, x), G_io,τ (:857-887, Theorem 4.1, bậc 3);
  - giao thức đầy đủ 18 bước (:1039-1061): sum-check #1 (µ1 = log m, ℓ1 = 3,
    T1 = 0), Prover gửi v_A, v_B, v_C, Verifier kiểm e_x = (v_A v_B - v_C)·eq(r_x, τ),
    rút r_A, r_B, r_C, sum-check #2 (ℓ2 = 2) cho
        M_{r_x}(y) = (r_A·A~(r_x,y) + r_B·B~(r_x,y) + r_C·C~(r_x,y))·Z~(y)   (:998-1005),
    cuối cùng mở cam kết tại r_y và kiểm e_y = (r_A v_1 + r_B v_2 + r_C v_3)·v_Z.

Khác biệt so với bài gốc (phải nêu trong sách):
  1. z = (1, x, w) theo quy ước snarklib, không phải Z = (io, 1, w) (:804).
     Spartan tách Z~(r_y) = (1 - r_y[0])·w~ + r_y[0]·(io,1)~ (:1027-1030); ở đây ta cam
     kết bảng W = z với các vị trí công khai (hằng 1 và x) đặt bằng 0, và Verifier
     tự cộng phần công khai: z~(r_y) = W~(r_y) + Σ_{i công khai} z_i·eq(i, r_y).
  2. PCS là snarklib.ligero (mã RS + Merkle), không phải Hyrax/KZG.
  3. Verifier tự tính A~(r_x, r_y) trực tiếp từ ma trận thưa trong O(nnz) —
     KHÔNG dùng SPARK (:1223 trở đi), nên Verifier không succinct theo kích thước mạch.
  4. Không có zero-knowledge (không che sum-check).
Ma trận được đệm: m (số ràng buộc) và n (số biến) lên lũy thừa của 2 và đệm
thêm để m = n (= 2^s) cho gọn (bài gốc cũng giả thiết m = n sau khi đệm).

API
    inst = pad_r1cs(r1cs)            -> dict A, B, C dạng thưa {(i, j): hệ số}, s, n_public
    proof = prove(inst, z, transcript, pcs=ligero, t=8)
    ok = verify(inst, x, proof, transcript, pcs=ligero, t=8)
    prove(..., pcs=None) chạy phần IOP (không cam kết), dùng cho đếm phép nhân.
"""

from .mle import eq_table, eq, evaluate, num_vars
from . import sumcheck as sc


def _pow2(n):
    s = 1
    while s < n:
        s *= 2
    return s


def pad_r1cs(r1cs):
    """Đệm R1CS (snarklib.r1cs) thành ma trận vuông 2^s × 2^s, lưu thưa."""
    size = max(_pow2(r1cs.m), _pow2(r1cs.n))
    mats = []
    for M in (r1cs.A, r1cs.B, r1cs.C):
        sp = {}
        for i, row in enumerate(M):
            for j, a in enumerate(row):
                if int(a) % r1cs.F.p:
                    sp[(i, j)] = int(a) % r1cs.F.p
        mats.append(sp)
    return {"A": mats[0], "B": mats[1], "C": mats[2], "size": size,
            "s": num_vars([0] * size), "n_public": r1cs.n_public,
            "m": r1cs.m, "n": r1cs.n}


def pad_z(inst, z):
    return list(z) + [0] * (inst["size"] - len(z))


def matvec_sparse(F, M, z, size):
    out = [F(0)] * size
    for (i, j), a in M.items():
        out[i] = out[i] + z[j] * a
    return out


def row_bind(F, M, ex, size):
    """Bảng y -> M~(r_x, y) = Σ_i eq(r_x, i)·M[i][y] (ex = eq_table(r_x)), thưa."""
    out = [F(0)] * size
    for (i, j), a in M.items():
        out[j] = out[j] + ex[i] * a
    return out


def matrix_eval(F, M, ex, ey):
    """A~(r_x, r_y) = Σ_{(i,j)} A_ij·eq(r_x,i)·eq(r_y,j)  — O(nnz), không SPARK."""
    return sum((ex[i] * ey[j] * a for (i, j), a in M.items()), F(0))


def witness_table(inst, z):
    """W = z với các vị trí công khai (0..n_public) đặt 0 — bảng được cam kết."""
    k = 1 + inst["n_public"]
    return [0] * k + list(z[k:])


def _g1(u):
    """eq·(a·b - c), bậc 3."""
    return u[0] * (u[1] * u[2] - u[3])


def _g2(u):
    return u[0] * u[1]


def prove(inst, z, transcript, pcs=None, t=8, F=None, label=b"spartan", phases=None):
    """z: véc-tơ (1, x, w) chưa đệm hoặc đã đệm. Trả về proof (dict).

    phases (dict, tùy chọn) nhận số phép nhân từng pha khi F là trường đếm
    (sumcheck.counting_field): matvec, eq_tau, sc1, row_bind, sc2, eval_w.
    """
    if F is None:
        F = pcs.F
    mark = [sc.mul_count(F)]

    def tick(name):
        if phases is not None:
            now = sc.mul_count(F)
            phases[name] = phases.get(name, 0) + now - mark[0]
            mark[0] = now
    size, s = inst["size"], inst["s"]
    z = [F(v) for v in pad_z(inst, z)]
    W = witness_table(inst, z)
    proof = {}
    if pcs is not None:
        root, st = pcs.commit([int(v) for v in W])
        proof["root"] = root
        transcript.absorb(label + b"/root", root)
    transcript.absorb(label + b"/x", [int(v) for v in z[1:1 + inst["n_public"]]])
    tau = [sc._challenge(F, transcript, label + b"/tau%d" % i) for i in range(s)]
    # sum-check #1: Σ_x eq(τ,x)·(Az~(x)·Bz~(x) - Cz~(x)) = 0
    Az = matvec_sparse(F, inst["A"], z, size)
    Bz = matvec_sparse(F, inst["B"], z, size)
    Cz = matvec_sparse(F, inst["C"], z, size)
    tick("matvec")
    E = eq_table(F, tau)
    tick("eq_tau")
    p1, rx, fin1 = sc.prove(F, [E, Az, Bz, Cz], _g1, 3, transcript, claim=F(0), label=label + b"/sc1")
    tick("sc1")
    vA, vB, vC = fin1[1], fin1[2], fin1[3]
    transcript.absorb(label + b"/vABC", [int(vA), int(vB), int(vC)])
    rA, rB, rC = [sc._challenge(F, transcript, label + b"/r" + c) for c in (b"A", b"B", b"C")]
    # sum-check #2: Σ_y (rA·A~(rx,y) + rB·B~(rx,y) + rC·C~(rx,y))·z~(y) = rA·vA + rB·vB + rC·vC
    ex = eq_table(F, rx)
    tA, tB, tC = (row_bind(F, inst[k], ex, size) for k in ("A", "B", "C"))
    Mrx = [rA * a + rB * b + rC * c for a, b, c in zip(tA, tB, tC)]
    T2 = rA * vA + rB * vB + rC * vC
    tick("row_bind")
    p2, ry, fin2 = sc.prove(F, [Mrx, z], _g2, 2, transcript, claim=T2, label=label + b"/sc2")
    tick("sc2")
    vW = evaluate(F, [F(v) for v in W], ry)      # W~(r_y), mở qua PCS
    tick("eval_w")
    transcript.absorb(label + b"/vW", int(vW))
    proof.update({"sc1": p1, "vABC": (vA, vB, vC), "sc2": p2, "vW": vW})
    if pcs is not None:
        proof["pcs"] = pcs.open(st, [int(v) for v in ry], transcript, t=t)
    return proof


def verify(inst, x, proof, transcript, pcs, t=8, label=b"spartan", info=None):
    """info (dict, tùy chọn) nhận r_x, r_y và các giá trị trung gian để in minh họa."""
    F = pcs.F
    size, s, k = inst["size"], inst["s"], inst["n_public"]
    if len(x) != k:
        return False
    if "root" not in proof:
        return False
    transcript.absorb(label + b"/root", proof["root"])
    transcript.absorb(label + b"/x", [int(v) for v in x])
    tau = [sc._challenge(F, transcript, label + b"/tau%d" % i) for i in range(s)]
    ok, rx, ex_ = sc.verify(F, s, 3, proof["sc1"], transcript, label + b"/sc1", claim=0)
    if not ok:
        return False
    vA, vB, vC = (F(v) for v in proof["vABC"])
    if (vA * vB - vC) * eq(F, rx, tau) != ex_:
        return False
    transcript.absorb(label + b"/vABC", [int(vA), int(vB), int(vC)])
    rA, rB, rC = [sc._challenge(F, transcript, label + b"/r" + c) for c in (b"A", b"B", b"C")]
    T2 = rA * vA + rB * vB + rC * vC
    ok, ry, ey = sc.verify(F, s, 2, proof["sc2"], transcript, label + b"/sc2", claim=T2)
    if not ok:
        return False
    vW = F(proof["vW"])
    transcript.absorb(label + b"/vW", int(vW))
    if F(proof["pcs"]["v"]) != vW:
        return False
    if not pcs.verify(proof["root"], [int(v) for v in ry], vW, proof["pcs"], transcript, t=t):
        return False
    ey_tab = eq_table(F, ry)
    pub = [F(1)] + [F(v) for v in x]
    vZ = vW + sum((pub[i] * ey_tab[i] for i in range(1 + k)), F(0))
    exr = eq_table(F, rx)
    v1, v2, v3 = (matrix_eval(F, inst[m], exr, ey_tab) for m in ("A", "B", "C"))
    if info is not None:
        info.update({"tau": tau, "rx": rx, "ry": ry, "vABC": (vA, vB, vC),
                     "rABC": (rA, rB, rC), "vZ": vZ, "v123": (v1, v2, v3), "ex": ex_, "ey": ey})
    return (rA * v1 + rB * v2 + rC * v3) * vZ == ey


def proof_size(inst, pcs, t=8):
    """Đếm phần tử trường và giá trị băm trong proof (sum-check gửi đủ d+1 giá trị/vòng)."""
    s = inst["s"]
    sz = pcs.sizes(inst["size"], t=t)
    sc_field = s * 4 + s * 3          # khẳng định 0 và T2 Verifier tự biết, không gửi
    return {"sumcheck1": s * 4, "sumcheck2": s * 3, "vABC": 3, "vW": 1,
            "sumcheck_optimized": s * 3 + s * 2,     # bỏ s_j(1) (Verifier tự suy)
            "pcs_field": sz["proof_field"], "pcs_hash": sz["proof_hash"], "commit_hash": 1,
            "rounds": 2 * s, "queries": sz["t"], "m_ligero": sz["m"], "k_ligero": sz["k"],
            "field_total": sc_field + 3 + 1 + sz["proof_field"],
            "hash_total": sz["proof_hash"] + 1}


def prover_mul_count(F_count, inst, z, transcript):
    """(tổng, theo pha) số phép nhân trường của phần IOP (không PCS) trên trường đếm."""
    F_count.counter["mul"] = 0
    phases = {}
    prove(inst, z, transcript, pcs=None, F=F_count, phases=phases)
    return F_count.counter["mul"], phases
