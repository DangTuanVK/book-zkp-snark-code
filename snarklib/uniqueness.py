"""Kiểm tính duy nhất của witness cho R1CS nhỏ (Chương 12: an toàn mạch).

Ba công cụ, từ chắc chắn nhất tới rẻ nhất:

1. Vét cạn có lan truyền (enumerate_witnesses): liệt kê MỌI witness w thỏa
   (Az)∘(Bz) = Cz với z = (1, x, w), x cố định, trên toàn bộ F_p^{|w|}.
   Không cần biết mạch được viết thế nào: mỗi ràng buộc chỉ còn một ẩn u trở
   thành phương trình bậc <= 2 theo u, nên tập nghiệm của u hữu hạn và tính
   được; ràng buộc bit u·u = u cho ngay {0, 1}. Khi không lan truyền tiếp được
   nữa, rẽ nhánh trên biến có miền nhỏ nhất. Kết quả là đầy đủ (không bỏ sót)
   miễn là số nút tìm kiếm không vượt ngân sách.

2. SMT (find_witness_z3, second_witness_z3, all_witnesses_z3): mã hóa mỗi
   biến là bit-vector 0 <= v < p và mỗi ràng buộc là
   (A·z mod p)·(B·z mod p) mod p = C·z mod p; hỏi z3 "có witness khác w0
   không?".
   Module import z3 có điều kiện (HAVE_Z3); phần 1 và 3 chạy bằng Python thuần.

3. Fuzz (fuzz_witness): đột biến một witness hợp lệ và đếm số đột biến vẫn
   thỏa hệ. Rẻ nhưng không đầy đủ: fuzz không thấy lỗi không có nghĩa là mạch
   an toàn.
"""
import random

try:                                   # z3 là tùy chọn
    import z3
    HAVE_Z3 = True
except ImportError:                    # pragma: no cover
    z3 = None
    HAVE_Z3 = False


# ---------------------------------------------------------------------------
# Biểu diễn thưa của R1CS
# ---------------------------------------------------------------------------
def sparse_rows(r1cs):
    """Trả về danh sách ràng buộc [(A_i, B_i, C_i)], mỗi LC là dict {cột: hệ số}."""
    p = r1cs.F.p

    def sp(row):
        return {j: int(c) % p for j, c in enumerate(row) if int(c) % p}

    return [(sp(a), sp(b), sp(c)) for a, b, c in zip(r1cs.A, r1cs.B, r1cs.C)]


def _split(lc, z, p):
    """LC tại phép gán từng phần z (None = chưa biết) -> (hằng, {ẩn: hệ số})."""
    k, unk = 0, {}
    for j, c in lc.items():
        if z[j] is None:
            unk[j] = c
        else:
            k = (k + c * z[j]) % p
    return k, unk


def _roots(p, c2, c1, c0, dom):
    """Nghiệm của c2·u^2 + c1·u + c0 = 0 trong dom (None = cả F_p).

    Trả về None nếu phương trình đồng nhất 0 (không thu hẹp được miền)."""
    c2, c1, c0 = c2 % p, c1 % p, c0 % p
    if c2 == 0 and c1 == 0:
        return None if c0 == 0 else set()
    if c2 == 0:
        u = (-c0 * pow(c1, p - 2, p)) % p
        return {u} if dom is None or u in dom else set()
    cand = range(p) if dom is None else dom
    return {u for u in cand if (c2 * u * u + c1 * u + c0) % p == 0}


def _propagate(rows, z, dom, p):
    """Lan truyền tới điểm bất động. Trả về False nếu mâu thuẫn."""
    changed = True
    while changed:
        changed = False
        for A, B, C in rows:
            a0, au = _split(A, z, p)
            b0, bu = _split(B, z, p)
            c0, cu = _split(C, z, p)
            unk = set(au) | set(bu) | set(cu)
            if not unk:
                if (a0 * b0 - c0) % p:
                    return False
                continue
            if len(unk) != 1:
                continue
            (u,) = unk
            a1, b1, c1 = au.get(u, 0), bu.get(u, 0), cu.get(u, 0)
            # (a0 + a1 u)(b0 + b1 u) - (c0 + c1 u) = 0
            R = _roots(p, a1 * b1, a0 * b1 + a1 * b0 - c1, a0 * b0 - c0, dom.get(u))
            if R is None:
                continue
            if not R:
                return False
            if len(R) == 1:
                z[u] = next(iter(R))
                dom.pop(u, None)
                changed = True
            elif dom.get(u) is None or len(R) < len(dom[u]):
                dom[u] = R
                changed = True
    return True


def enumerate_witnesses(r1cs, x, max_nodes=10 ** 6, stats=None):
    """Liệt kê mọi w (list int) thỏa R1CS với đầu vào công khai x.

    stats (dict, tùy chọn) nhận: nodes (số nút tìm kiếm), branch_vars (các
    biến từng bị rẽ nhánh), complete (True nếu không chạm ngân sách).
    """
    p = r1cs.F.p
    rows = sparse_rows(r1cs)
    n, ell = r1cs.n, r1cs.n_public
    z0 = [1] + [int(v) % p for v in x] + [None] * (n - 1 - ell)
    out, nodes, branched = [], [0], set()

    def rec(z, dom):
        nodes[0] += 1
        if nodes[0] > max_nodes:
            raise RuntimeError("vượt ngân sách vét cạn")
        if not _propagate(rows, z, dom, p):
            return
        free = [j for j in range(n) if z[j] is None]
        if not free:
            out.append([z[j] for j in range(1 + ell, n)])
            return
        # rẽ nhánh trên biến có miền nhỏ nhất (miền None = cả F_p)
        j = min(free, key=lambda v: len(dom[v]) if dom.get(v) is not None else p)
        branched.add(j)
        vals = sorted(dom[j]) if dom.get(j) is not None else range(p)
        for v in vals:
            z2, d2 = list(z), {k: set(s) for k, s in dom.items()}
            z2[j] = v
            d2.pop(j, None)
            rec(z2, d2)

    rec(list(z0), {})
    out.sort()
    if stats is not None:
        stats.update(nodes=nodes[0], branch_vars=sorted(branched), complete=True)
    return out


def brute_force_inputs(circ, fixed, ranges, r1cs=None):
    """Vét cạn ở mức mạch: thử mọi tổ hợp giá trị của biến bí mật trong ranges.

    circ: snarklib.dsl.Circuit; fixed: {tên public: giá trị};
    ranges: {tên private: iterable}. Witness còn lại sinh bằng gợi ý (hint).
    Chỉ đầy đủ khi ranges phủ mọi giá trị khả dĩ (ví dụ [0, 31] khi biến đã
    được phân tích 5 bit). Trả về danh sách dict các tổ hợp thỏa.
    """
    r1cs = r1cs or circ.compile_r1cs()
    names = sorted(ranges)
    out = []

    def rec(i, cur):
        if i == len(names):
            inp = dict(fixed, **cur)
            x, w = circ.r1cs_assignment(inp)
            if r1cs.is_satisfied(x, w):
                out.append(dict(cur))
            return
        for v in ranges[names[i]]:
            cur[names[i]] = v
            rec(i + 1, cur)

    rec(0, {})
    return out


def witness_names(circ):
    """Tên các phần tử của w theo thứ tự z = (1, x, w) của Circuit.compile_r1cs."""
    order = circ.z_order()
    ell = sum(1 for k in circ.kind if k == "public")
    return [circ.names[v] for v in order[1 + ell:]]


# ---------------------------------------------------------------------------
# SMT với z3 (tùy chọn)
# ---------------------------------------------------------------------------
Z3_WIDTH = 32   # độ rộng bit-vector: đủ chứa Σ|c|·(p-1) và tích hai số < p


def _z3_model(r1cs, x, extra=None, width=Z3_WIDTH):
    """Tạo (solver, danh sách biến z3 của w). extra(s, w) thêm ràng buộc riêng.

    Mã hóa: mỗi biến là bit-vector width bit với 0 <= v < p; mỗi LC được rút gọn
    (URem p) trước khi nhân nên mọi giá trị trung gian nằm gọn trong width bit
    (kiểm bằng assert). Bản số nguyên Int + mod (số học phi tuyến) z3 trả về
    "unknown" ngay trên mạch 13 ràng buộc, nên ta dùng bit-vector (quyết định
    được bằng bit-blasting).
    """
    if not HAVE_Z3:
        raise RuntimeError("không có z3: cài z3-solver hoặc dùng enumerate_witnesses")
    p = r1cs.F.p
    ell = r1cs.n_public
    rows = sparse_rows(r1cs)
    worst = max((sum(L.values()) for row in rows for L in row), default=0) * (p - 1)
    assert worst < 2 ** width and (p - 1) ** 2 < 2 ** width, "tăng width"
    BV = lambda v: z3.BitVecVal(v, width)
    w = [z3.BitVec("w%d" % i, width) for i in range(r1cs.n - 1 - ell)]
    z = [BV(1)] + [BV(int(v) % p) for v in x] + w
    P = BV(p)
    s = z3.Solver()
    for v in w:
        s.add(z3.ULT(v, P))

    def lc(L):
        if not L:
            return BV(0)
        return z3.URem(z3.Sum([BV(c) * z[j] for j, c in L.items()]), P)

    for A, B, C in rows:
        s.add(z3.URem(lc(A) * lc(B), P) == lc(C))
    if extra:
        extra(s, w)
    return s, w


def find_witness_z3(r1cs, x, exclude=(), extra=None, timeout_ms=60000):
    """Một witness thỏa R1CS khác mọi witness trong exclude, hoặc None (unsat)."""
    s, w = _z3_model(r1cs, x, extra)
    s.set("timeout", timeout_ms)
    for old in exclude:
        s.add(z3.Or([v != int(o) for v, o in zip(w, old)]))
    res = s.check()
    if res == z3.unsat:
        return None
    if res != z3.sat:
        raise RuntimeError("z3 trả về %s" % res)
    m = s.model()
    out = [m.eval(v, model_completion=True).as_long() for v in w]
    if not r1cs.is_satisfied(x, out):          # kiểm chéo mô hình của z3
        raise AssertionError("mô hình z3 không thỏa R1CS")
    return out


def second_witness_z3(r1cs, x, w0, **kw):
    """Câu hỏi kiểu Picus: có witness w != w0 cùng đầu vào x không? (None = không)."""
    return find_witness_z3(r1cs, x, exclude=[[int(v) for v in w0]], **kw)


def all_witnesses_z3(r1cs, x, limit=64, **kw):
    """Lặp find_witness_z3 với mệnh đề chặn tới khi unsat; trả về danh sách đã sắp."""
    found = []
    while len(found) < limit:
        w = find_witness_z3(r1cs, x, exclude=found, **kw)
        if w is None:
            break
        found.append(w)
    return sorted(found)


# ---------------------------------------------------------------------------
# Fuzz witness
# ---------------------------------------------------------------------------
MUTATORS = ("random_one", "plus_minus_one", "swap_two", "flip_bit", "negate")


def mutate(w, kind, rng, p):
    """Một đột biến của w (list int). Trả về list mới."""
    w = list(w)
    i = rng.randrange(len(w))
    if kind == "random_one":            # thay một ô bằng giá trị ngẫu nhiên
        w[i] = rng.randrange(p)
    elif kind == "plus_minus_one":      # ±1 một ô
        w[i] = (w[i] + rng.choice((1, -1))) % p
    elif kind == "swap_two":            # tráo hai ô
        j = rng.randrange(len(w))
        w[i], w[j] = w[j], w[i]
    elif kind == "flip_bit":            # lật bit thấp nhất của một ô
        w[i] = w[i] ^ 1 if (w[i] ^ 1) < p else w[i]
    elif kind == "negate":              # đổi dấu một ô: v -> p - v
        w[i] = (-w[i]) % p
    else:
        raise ValueError(kind)
    return w


def fuzz_witness(r1cs, x, w, trials=10000, rng=None, kinds=MUTATORS):
    """Đột biến witness hợp lệ w ở mức R1CS; đếm số đột biến vẫn thỏa.

    Trả về {kind: (số lần thử, số đột biến thỏa và KHÁC w, tập witness mới)}.
    """
    rng = rng or random.Random(0)
    p = r1cs.F.p
    w = [int(v) % p for v in w]
    if not r1cs.is_satisfied(x, w):
        raise ValueError("witness gốc không thỏa")
    res = {}
    for kind in kinds:
        hit, new = 0, set()
        for _ in range(trials):
            m = mutate(w, kind, rng, p)
            if m != w and r1cs.is_satisfied(x, m):
                hit += 1
                new.add(tuple(m))
        res[kind] = (trials, hit, new)
    return res


def fuzz_inputs(circ, fixed, base, trials=10000, rng=None, spread=None, r1cs=None):
    """Fuzz ở mức mạch: đột biến các biến bí mật, sinh lại witness bằng gợi ý.

    base: {tên private: giá trị hợp lệ}; spread: {tên: iterable các giá trị thử}
    (mặc định cả F_p). Trả về (số lần thử, {tổ hợp thỏa khác base: số lần gặp}).
    """
    rng = rng or random.Random(0)
    r1cs = r1cs or circ.compile_r1cs()
    p = circ.F.p
    names = sorted(base)
    spread = spread or {k: range(p) for k in names}
    pools = {k: list(spread[k]) for k in names}
    seen = {}
    for _ in range(trials):
        cur = {k: rng.choice(pools[k]) for k in names}
        if cur == base:
            continue
        x, w = circ.r1cs_assignment(dict(fixed, **cur))
        if r1cs.is_satisfied(x, w):
            key = tuple(cur[k] for k in names)
            seen[key] = seen.get(key, 0) + 1
    return trials, seen
