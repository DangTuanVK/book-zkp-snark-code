"""Chương 5 (mã phụ) -- Nhúng Plonkish và AIR vào CCS trên mạch N = 391.

Số liệu cho Mục 5.1 (CCS):
  (1) Plonkish 20 hàng của snarklib.dsl -> CCS: mỗi ô dây ứng với một cột của z
      (ràng buộc sao chép được hấp thụ vì hai ô nối nhau dùng chung một cột),
      selector gập vào hệ số của ma trận (mỗi đơn thức của cổng chứa đúng một
      selector). Bản "tách" giữ mỗi đơn thức một đa tập; bản "gộp" dồn các đơn
      thức bậc 1 vào một ma trận.
  (2) AIR 6 hàng x 5 cột cho cùng mệnh đề (bit cao trước, tích lũy Horner),
      in trace của (17, 23) và kiểm bằng CCS.
  (3) Cổng bậc cao: y = x^5 là một hàng CCS với đa tập {0,0,0,0,0}, so với 3
      ràng buộc R1CS.
Quy ước z = (1, x, w) như snarklib.ccs (bài CCS gốc dùng (w, 1, x)).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.field import GF
from snarklib.ccs import CCS
from snarklib.dsl import mach_phan_tich

F = GF(449)
p = F.p
show = lambda v: int(v) - p if int(v) > p // 2 else int(v)


def nnz(M):
    return sum(1 for row in M for v in row if int(v) % p)


# ---------------------------------------------------------------- (1) Plonkish -> CCS
def plonkish_to_ccs(c, gop):
    """CCS từ mạch Plonkish 3 cột đã biên dịch bởi Circuit c."""
    pk = c.compile_plonkish()
    cells, temps, n_rows = c._plonk
    n_pub = pk.n_public
    # cột của z: 0 = hằng 1, 1..n_pub = x, sau đó các biến bí mật và dây tạm
    col = {}
    pubs = [v for v, k in enumerate(c.kind) if k == "public"]
    for i, v in enumerate(pubs):
        col[("v", v)] = 1 + i
    nxt = 1 + n_pub
    for ws in cells:
        for w in ws:
            if w is not None and w not in col:
                col[w] = nxt
                nxt += 1
    n = nxt
    rows = list(range(n_pub, n_rows))   # hàng công khai trở thành tầm thường
    m = len(rows)

    def zero():
        return [[0] * n for _ in range(m)]

    sel = pk.sel
    if gop:
        M = [zero(), zero(), zero()]    # M0 = q_M * a, M1 = b, M2 = tuyến tính
        S, cst = [[0, 1], [2]], [1, 1]
    else:
        M = [zero() for _ in range(6)]  # q_M*a, b, q_L*a, q_R*b, q_O*c, q_C*1
        S, cst = [[0, 1], [2], [3], [4], [5]], [1, 1, 1, 1, 1]
    for r, j in enumerate(rows):
        wa, wb, wc = cells[j]
        qL, qR, qO, qM, qC = (int(sel[s][j]) for s in ("qL", "qR", "qO", "qM", "qC"))
        if qM:
            M[0][r][col[wa]] = qM
            M[1][r][col[wb]] = 1
        lin = [(2, wa, qL), (3, wb, qR), (4, wc, qO)]
        for k, w, q in lin:
            if q:
                tgt = M[2] if gop else M[k]
                tgt[r][col[w]] = (tgt[r][col[w]] + q) % p
        if qC:
            tgt = M[2] if gop else M[5]
            tgt[r][0] = (tgt[r][0] + qC) % p
    ccs = CCS(F, M, S, cst, n_pub)
    return ccs, col, rows, pk


def ccs_assignment(c, col, inputs):
    cells, temps, _ = c._plonk
    env = c.solve(inputs)
    n = max(col.values()) + 1
    z = [None] * n
    z[0] = F(1)
    for w, j in col.items():
        z[j] = env[w[1]] if w[0] == "v" else c.eval(temps[w[1]], env)
    n_pub = sum(1 for k in c.kind if k == "public")
    return z[1:1 + n_pub], z[1 + n_pub:]


c = mach_phan_tich()[0]
print("== (1) Plonkish -> CCS, mạch N = 391 ==")
for gop in (False, True):
    ccs, col, rows, pk = plonkish_to_ccs(c, gop)
    x, w = ccs_assignment(c, col, {"N": 391, "a": 17, "b": 23})
    ok = ccs.is_satisfied(x, w)
    x2, w2 = ccs_assignment(c, col, {"N": 391, "a": 17, "b": 24})
    bad = ccs.is_satisfied(x2, w2)
    print("%s: m = %d (bỏ %d hàng công khai), n = %d (1 + %d công khai + %d bí mật/tạm),"
          " t = %d, q = %d, d = %d, c = %s"
          % ("gộp " if gop else "tách", ccs.m, pk.n - ccs.m, ccs.n, ccs.n_public,
             ccs.n - 1 - ccs.n_public, len(ccs.M), len(ccs.S), ccs.degree(), [show(v) for v in ccs.c]))
    print("      nnz(M_j) = %s; (17, 23) thỏa: %s; (17, 24) thỏa: %s"
          % ([nnz(Mj) for Mj in ccs.M], ok, bad))
    assert ok and not bad
print("      %d ràng buộc sao chép được hấp thụ: %d ô dây -> %d cột witness"
      % (len(pk.copies), sum(1 for ws in c._plonk[0][1:pk.n] for wv in ws if wv is not None),
         ccs.n - 1 - ccs.n_public))
# hàng 19 (cổng nhân a·b = N) trong bản gộp
ccs, col, rows, pk = plonkish_to_ccs(c, True)
r = rows.index(19)
names = {j: (c.names[w[1]] if w[0] == "v" else "t%d" % w[1]) for w, j in col.items()}
names[0] = "1"
desc = lambda Mj: " + ".join("%s·%s" % (show(v), names[k]) for k, v in enumerate(Mj[r]) if v) or "0"
print("      hàng cổng nhân (hàng Plonkish 19): (M0 z)=%s, (M1 z)=%s, (M2 z)=%s"
      % (desc(ccs.M[0]), desc(ccs.M[1]), desc(ccs.M[2])))

# ---------------------------------------------------------------- (2) AIR
print("== (2) AIR 6 hàng x 5 cột (bit cao trước) ==")
K = 5


def trace(a, b):
    ab = [0] + [(a >> (K - 1 - i)) & 1 for i in range(K)]
    bb = [0] + [(b >> (K - 1 - i)) & 1 for i in range(K)]
    A, B = [0], [0]
    for i in range(1, K + 1):
        A.append((2 * A[-1] + ab[i]) % p)
        B.append((2 * B[-1] + bb[i]) % p)
    P = [(x * y) % p for x, y in zip(A, B)]
    return [ab, A, bb, B, P]       # cột: alpha, A, beta, B, P


def air_ok(T, N):
    al, A, be, B, P = T
    T_ = len(A)
    row = all(v * (v - 1) % p == 0 for v in al + be) and all((P[i] - A[i] * B[i]) % p == 0 for i in range(T_))
    tr = all((A[i + 1] - 2 * A[i] - al[i + 1]) % p == 0 and (B[i + 1] - 2 * B[i] - be[i + 1]) % p == 0
             for i in range(T_ - 1))
    bd = al[0] == A[0] == be[0] == B[0] == 0 and P[T_ - 1] % p == N % p
    return row and tr and bd


for (a, b) in ((17, 23), (17, 24)):
    T = trace(a, b)
    print("(a, b) = (%d, %d)" % (a, b))
    for nm, colv in zip(("alpha", "A", "beta", "B", "P"), T):
        print("   %-5s = %s" % (nm, colv))
    print("   thỏa AIR (N = 391):", air_ok(T, 391))
assert air_ok(trace(17, 23), 391) and not air_ok(trace(17, 24), 391)


def air_to_ccs(N_rows=K + 1):
    """CCS cho AIR: z = (1, N, các ô trace theo hàng); t = 3, S = ({0,1},{2})."""
    w_ = 5
    n = 2 + w_ * N_rows
    cell = lambda i, j: 2 + w_ * i + j          # j: 0 alpha, 1 A, 2 beta, 3 B, 4 P
    M0, M1, M2 = [], [], []

    def add(l=None, r=None, lin=()):
        a0, a1, a2 = [0] * n, [0] * n, [0] * n
        if l is not None:
            a0[l[0]] = l[1]
            a1[r] = 1
        for k, v in lin:
            a2[k] = (a2[k] + v) % p
        M0.append(a0); M1.append(a1); M2.append(a2)

    for i in range(N_rows):                     # ràng buộc hàng (mọi hàng)
        add((cell(i, 0), 1), cell(i, 0), [(cell(i, 0), -1)])        # alpha^2 - alpha
        add((cell(i, 2), 1), cell(i, 2), [(cell(i, 2), -1)])        # beta^2 - beta
        add((cell(i, 1), -1), cell(i, 3), [(cell(i, 4), 1)])        # P - A·B
    for i in range(N_rows - 1):                 # ràng buộc chuyển (quay +1)
        add(lin=[(cell(i + 1, 1), 1), (cell(i, 1), -2), (cell(i + 1, 0), -1)])
        add(lin=[(cell(i + 1, 3), 1), (cell(i, 3), -2), (cell(i + 1, 2), -1)])
    for j in range(4):                          # biên: hàng 0 bằng 0
        add(lin=[(cell(0, j), 1)])
    add(lin=[(cell(N_rows - 1, 4), 1), (1, -1)])  # biên: P cuối = N
    return CCS(F, [M0, M1, M2], [[0, 1], [2]], [1, 1], 1)


ccs_air = air_to_ccs()
T = trace(17, 23)
w_air = [T[j][i] for i in range(K + 1) for j in range(5)]
ok = ccs_air.is_satisfied([391], w_air)
print("AIR -> CCS: m = %d (%d hàng + %d chuyển + %d biên), n = %d, t = %d, q = %d, d = %d; (17, 23) thỏa: %s"
      % (ccs_air.m, 3 * (K + 1), 2 * K, 5, ccs_air.n, len(ccs_air.M), len(ccs_air.S), ccs_air.degree(), ok))
assert ok

# ---------------------------------------------------------------- (3) cổng bậc 5
print("== (3) Cổng bậc cao y = x^5 ==")
Mx = [[0, 0, 1]]         # z = (1, y, x): chọn x
My = [[0, 1, 0]]         # chọn y
ccs5 = CCS(F, [Mx, My], [[0, 0, 0, 0, 0], [1]], [1, -1], 1)
print("CCS: m = 1, t = 2, S = %s, c = (1, -1), d = %d" % (ccs5.S, ccs5.degree()))
print("x = 3, y = %d: thỏa %s; y = 244: thỏa %s"
      % (int(F(3) ** 5), ccs5.is_satisfied([int(F(3) ** 5)], [3]), ccs5.is_satisfied([244], [3])))
print("R1CS cùng cổng: x2 = x·x, x4 = x2·x2, y = x4·x -> 3 ràng buộc, 2 biến trung gian")
