"""DSL tí hon mô tả mạch số học và biên dịch sang R1CS, Plonkish, CCS.

Mô hình:
    c = Circuit()                 # mặc định trên F_449
    N = c.public("N")             # biến công khai (nằm trong x)
    a = c.private("a")            # biến bí mật (nằm trong w)
    bits = c.to_bits(a, 5)        # 5 biến bit + ràng buộc bit + ràng buộc tổ hợp
    c.assert_mul(a, b, N)         # a·b = N
    r1cs = c.compile_r1cs()       # R1CS của snarklib.r1cs, z = (1, x, w)
    x, w = c.r1cs_assignment({"N": 391, "a": 17, "b": 23})

Một biểu thức là tổ hợp tuyến tính (LC) của các biến: add/sub/scale/const
không sinh ràng buộc; chỉ mul, assert_* và to_bits sinh ràng buộc. Mỗi ràng
buộc lưu ở dạng R1CS tổng quát (LC_A)·(LC_B) = LC_C; trình biên dịch Plonkish
chẻ các LC dài thành chuỗi cổng cộng.

Biến trung gian (kết quả mul, bit, nghịch đảo) có "gợi ý" (hint) tính giá trị
từ các biến đã biết; solve() chạy các gợi ý theo thứ tự tạo để sinh witness.
Gợi ý không phải ràng buộc: witness sai (ví dụ a không vừa k bit) vẫn sinh ra
được, chỉ là không thỏa hệ ràng buộc.
"""

from .field import GF
from .r1cs import R1CS
from .plonkish import PlonkishCircuit
from .ccs import from_r1cs

ONE = 0  # chỉ số biến hằng 1


class LC:
    """Tổ hợp tuyến tính Σ c_v · v (v = 0 là hằng 1). Hệ số lưu dạng int mod p."""

    def __init__(self, circ, terms=None):
        self.circ = circ
        p = circ.F.p
        self.terms = {v: c % p for v, c in (terms or {}).items() if c % p}

    # --- số học tuyến tính (không sinh ràng buộc) ---
    def _lift(self, o):
        return o if isinstance(o, LC) else self.circ.const(o)

    def __add__(self, o):
        o = self._lift(o)
        t = dict(self.terms)
        for v, c in o.terms.items():
            t[v] = t.get(v, 0) + c
        return LC(self.circ, t)

    __radd__ = __add__

    def __neg__(self):
        return LC(self.circ, {v: -c for v, c in self.terms.items()})

    def __sub__(self, o):
        return self + (-self._lift(o))

    def __rsub__(self, o):
        return self._lift(o) - self

    def __mul__(self, k):
        """Chỉ nhân với hằng số; nhân hai biến dùng Circuit.mul."""
        if isinstance(k, LC):
            raise TypeError("nhân hai biểu thức: dùng Circuit.mul")
        k = int(k)
        return LC(self.circ, {v: c * k for v, c in self.terms.items()})

    __rmul__ = __mul__

    def constant(self):
        return self.terms.get(ONE, 0)

    def var_terms(self):
        return [(v, c) for v, c in sorted(self.terms.items()) if v != ONE]

    def is_constant(self):
        return not self.var_terms()

    def __repr__(self):
        names, p = self.circ.names, self.circ.F.p
        out = ""
        for v, c in sorted(self.terms.items(), key=lambda t: (t[0] == ONE, t[0])):
            c = c - p if c > p // 2 else c          # in hệ số âm cho dễ đọc
            mag = abs(c)
            term = str(mag) if v == ONE else (names[v] if mag == 1 else "%d·%s" % (mag, names[v]))
            out += ("-" if c < 0 else "") + term if not out else (" - " if c < 0 else " + ") + term
        return out or "0"


class Circuit:
    def __init__(self, F=None):
        self.F = F or GF(449)
        self.names = ["1"]          # tên biến, chỉ số 0 là hằng 1
        self.kind = ["one"]         # "one" | "public" | "private" | "aux"
        self.hints = {}             # biến aux -> hàm tính giá trị từ env
        self.constraints = []       # danh sách (LC_A, LC_B, LC_C, nhãn)

    # ------------------------------------------------------------ khai báo
    def _new_var(self, name, kind, hint=None):
        self.names.append(name)
        self.kind.append(kind)
        v = len(self.names) - 1
        if hint is not None:
            self.hints[v] = hint
        return LC(self, {v: 1})

    def public(self, name):
        return self._new_var(name, "public")

    def private(self, name):
        return self._new_var(name, "private")

    def const(self, c):
        return LC(self, {ONE: int(c)})

    def add(self, x, y):
        return x + y

    def sub(self, x, y):
        return x - y

    def scale(self, k, x):
        return x * k

    # ------------------------------------------------------------ ràng buộc
    def enforce(self, A, B, C, label=""):
        """Ràng buộc tổng quát A·B = C (A, B, C là LC hoặc hằng)."""
        lift = lambda e: e if isinstance(e, LC) else self.const(e)
        self.constraints.append((lift(A), lift(B), lift(C), label))

    def mul(self, x, y, name=None):
        """Trả về biến mới p với ràng buộc x·y = p (hằng số thì không sinh ràng buộc)."""
        if isinstance(x, LC) and x.is_constant():
            return y * x.constant()
        if isinstance(y, LC) and y.is_constant():
            return x * y.constant()
        name = name or "t%d" % len(self.names)
        p = self._new_var(name, "aux", hint=lambda env: self.eval(x, env) * self.eval(y, env))
        self.enforce(x, y, p, "mul")
        return p

    def assert_mul(self, x, y, z, label="mul"):
        """x·y = z (một ràng buộc, không sinh biến mới)."""
        self.enforce(x, y, z, label)

    def assert_eq(self, x, y, label="eq"):
        """x = y, biểu diễn (x - y)·1 = 0."""
        self.enforce(x - y, self.const(1), self.const(0), label)

    def assert_bool(self, x, label="bool"):
        """x ∈ {0, 1}, biểu diễn x·x = x."""
        self.enforce(x, x, x, label)

    def assert_nonzero(self, x, name=None, label="nonzero"):
        """x ≠ 0: thêm biến u với gợi ý u = x^{-1} và ràng buộc x·u = 1."""
        name = name or "inv%d" % len(self.names)

        def inv(env):
            v = self.eval(x, env)
            return v.inverse() if v != 0 else self.F(0)

        u = self._new_var(name, "aux", hint=inv)
        self.enforce(x, u, self.const(1), label)
        return u

    def to_bits(self, x, k, prefix=None):
        """Phân tích x thành k bit b_0..b_{k-1} (b_0 thấp nhất).

        Sinh k ràng buộc bit b_i·b_i = b_i và một ràng buộc tổ hợp
        (Σ 2^i b_i)·1 = x. Trả về danh sách k biểu thức bit.
        """
        if prefix is None:
            vt = x.var_terms()
            prefix = self.names[vt[0][0]] + "_" if len(x.terms) == 1 else "bit%d_" % len(self.names)
        bits = []
        for i in range(k):
            b = self._new_var("%s%d" % (prefix, i), "aux",
                              hint=lambda env, i=i: self.F((int(self.eval(x, env)) >> i) & 1))
            self.assert_bool(b, "bool")
            bits.append(b)
        total = self.const(0)
        for i, b in enumerate(bits):
            total = total + b * (1 << i)
        self.enforce(total, self.const(1), x, "recompose")
        return bits

    def assert_le(self, x, bound, k=None, prefix=None, label="le"):
        """x <= bound (bound là hằng số nguyên >= 0), thêm ở Chương 12.

        Cách làm: phân tích d = bound - x thành k bit (k = số bit của bound),
        tức k ràng buộc bit + 1 ràng buộc tổ hợp. Ràng buộc chỉ bảo đảm
        d ∈ [0, 2^k - 1], nghĩa là x ∈ {bound - d mod p}; muốn có đúng
        0 <= x <= bound thì x phải ĐÃ được chặn x ∈ [0, 2^kx - 1] (ví dụ bằng
        to_bits) và 2^kx + 2^k <= p + 1 để hai khoảng không quấn quanh p.
        Trả về danh sách k bit của d.
        """
        bound = int(bound)
        if bound < 0:
            raise ValueError("bound phải >= 0")
        k = k or max(1, bound.bit_length())
        if (1 << k) - 1 < bound:
            raise ValueError("k bit không chứa được bound")
        prefix = prefix or "le%d_" % len(self.names)
        d = self.const(bound) - x
        bits = []
        for i in range(k):
            b = self._new_var("%s%d" % (prefix, i), "aux",
                              hint=lambda env, i=i: self.F((int(self.eval(d, env)) >> i) & 1))
            self.assert_bool(b, "bool")
            bits.append(b)
        total = self.const(0)
        for i, b in enumerate(bits):
            total = total + b * (1 << i)
        self.enforce(total, self.const(1), d, label)
        return bits

    # ------------------------------------------------------------ witness
    def eval(self, lc, env):
        F = self.F
        return sum((F(c) * env[v] for v, c in lc.terms.items()), F(0))

    def solve(self, inputs):
        """inputs: {tên: giá trị} cho mọi biến public/private. Trả về env (list)."""
        F = self.F
        env = [None] * len(self.names)
        env[ONE] = F(1)
        for v, k in enumerate(self.kind):
            if k in ("public", "private"):
                env[v] = F(inputs[self.names[v]])
            elif k == "aux":
                env[v] = F(self.hints[v](env))
        return env

    def is_satisfied_env(self, env):
        return all(self.eval(A, env) * self.eval(B, env) == self.eval(C, env)
                   for A, B, C, _ in self.constraints)

    # ------------------------------------------------------------ R1CS
    def z_order(self):
        """Thứ tự biến trong z = (1, x, w): hằng 1, các biến public, còn lại."""
        pub = [v for v, k in enumerate(self.kind) if k == "public"]
        rest = [v for v, k in enumerate(self.kind) if k in ("private", "aux")]
        return [ONE] + pub + rest

    def compile_r1cs(self):
        order = self.z_order()
        pos = {v: j for j, v in enumerate(order)}
        n = len(order)

        def row(lc):
            r = [0] * n
            for v, c in lc.terms.items():
                r[pos[v]] = c
            return r

        A = [row(a) for a, _, _, _ in self.constraints]
        B = [row(b) for _, b, _, _ in self.constraints]
        C = [row(c) for _, _, c, _ in self.constraints]
        n_pub = sum(1 for k in self.kind if k == "public")
        return R1CS(self.F, A, B, C, n_public=n_pub)

    def r1cs_assignment(self, inputs):
        """Trả về (x, w) theo thứ tự z = (1, x, w)."""
        env = self.solve(inputs)
        order = self.z_order()
        n_pub = sum(1 for k in self.kind if k == "public")
        z = [env[v] for v in order]
        return z[1:1 + n_pub], z[1 + n_pub:]

    def public_values(self, inputs):
        return [self.F(inputs[self.names[v]]) for v, k in enumerate(self.kind) if k == "public"]

    # ------------------------------------------------------------ CCS
    def compile_ccs(self):
        """CCS bậc 2 nhúng từ R1CS: M = (A, B, C), S = ({0,1}, {2}), c = (1, -1)."""
        return from_r1cs(self.compile_r1cs())

    # ------------------------------------------------------------ Plonkish
    def compile_plonkish(self, pad_pow2=False):
        """Biên dịch sang mạch Plonkish 3 cột (snarklib.plonkish).

        Bố cục: hàng 0..ℓ-1 là cổng công khai (q_L = 1, a = x_i, PI_i = -x_i);
        sau đó mỗi ràng buộc thành một hoặc nhiều hàng. LC có >= 2 biến được gộp
        dần bằng cổng cộng s = c1·u + c2·v vào dây tạm. Cùng một biến xuất hiện
        ở nhiều ô thì nối bằng ràng buộc sao chép. Trả về PlonkishCircuit; dùng
        plonkish_assignment() để lấy ba cột a, b, c.
        """
        F = self.F
        p = F.p
        rows, cells = [], []      # cells[j] = (wa, wb, wc); dây = ("v", id) | ("t", id) | None
        temps = []                # temps[i] = LC (giá trị của dây tạm)

        def emit(sel, wa=None, wb=None, wc=None):
            rows.append({k: v % p for k, v in sel.items() if v % p})
            cells.append((wa, wb, wc))

        def wire_lc(w):
            return LC(self, {w[1]: 1}) if w[0] == "v" else temps[w[1]]

        def collapse(terms):
            """Gộp [(dây, hệ số), ...] (>= 2 phần tử) thành một dây tạm, trả về dây đó."""
            (w1, c1), (w2, c2) = terms[0], terms[1]
            for k in range(1, len(terms)):
                if k > 1:
                    w2, c2 = terms[k]
                t = ("t", len(temps))
                temps.append(wire_lc(w1) * c1 + wire_lc(w2) * c2)
                emit({"qL": c1, "qR": c2, "qO": -1}, w1, w2, t)
                w1, c1 = t, 1
            return w1

        def linear_zero(terms, k0):
            """Σ c·dây + k0 = 0."""
            terms = list(terms)
            while len(terms) > 3:
                s = collapse(terms[:2])
                terms = [(s, 1)] + terms[2:]
            ws = [w for w, _ in terms] + [None] * (3 - len(terms))
            cs = [c for _, c in terms] + [0] * (3 - len(terms))
            emit({"qL": cs[0], "qR": cs[1], "qO": cs[2], "qC": k0}, *ws)

        def single(lc):
            """LC -> (hệ số, dây, hằng) với tối đa một dây (gộp nếu cần)."""
            vt = [(("v", v), c) for v, c in lc.var_terms()]
            if not vt:
                return 0, None, lc.constant()
            if len(vt) == 1:
                return vt[0][1], vt[0][0], lc.constant()
            return 1, collapse(vt), lc.constant()

        pubs = [v for v, k in enumerate(self.kind) if k == "public"]
        for v in pubs:
            emit({"qL": 1}, ("v", v), None, None)

        for A, B, C, _ in self.constraints:
            if B.is_constant() or A.is_constant():
                if A.is_constant():
                    A, B = B, A
                lin = A * B.constant() - C          # A·b0 - C = 0
                vt = [(("v", v), c) for v, c in lin.var_terms()]
                linear_zero(vt, lin.constant())
                continue
            al, x, a0 = single(A)
            be, y, b0 = single(B)
            ga, zw, c0 = single(C)
            # (al·x + a0)(be·y + b0) - (ga·z + c0) = 0
            emit({"qM": al * be, "qL": al * b0, "qR": a0 * be, "qO": -ga,
                  "qC": a0 * b0 - c0}, x, y, zw)

        n_rows = len(rows)
        if pad_pow2:
            size = 1
            while size < n_rows:
                size *= 2
            for _ in range(size - n_rows):
                emit({}, None, None, None)

        # ràng buộc sao chép: nối mọi ô cùng một dây vào ô đầu tiên
        occ = {}
        for j, ws in enumerate(cells):
            for col, w in zip("abc", ws):
                if w is not None:
                    occ.setdefault(w, []).append((col, j))
        copies = []
        for w, lst in occ.items():
            for cell in lst[1:]:
                copies.append((lst[0], cell))

        self._plonk = (cells, temps, n_rows)
        circ = PlonkishCircuit(F, rows, copies, n_public=len(pubs))
        circ.n_used = n_rows
        circ.n_temps = len(temps)
        return circ

    def plonkish_assignment(self, inputs):
        """Ba cột (a, b, c) và x cho bản Plonkish đã biên dịch gần nhất."""
        cells, temps, _ = self._plonk
        env = self.solve(inputs)
        F = self.F

        def val(w):
            if w is None:
                return F(0)
            return env[w[1]] if w[0] == "v" else self.eval(temps[w[1]], env)

        cols = [[val(ws[k]) for ws in cells] for k in range(3)]
        return cols[0], cols[1], cols[2], self.public_values(inputs)

    # ------------------------------------------------------------ thống kê
    def stats(self):
        r = self.compile_r1cs()
        nnz = [sum(1 for row in M for v in row if v) for M in (r.A, r.B, r.C)]
        labels = {}
        for _, _, _, lab in self.constraints:
            labels[lab] = labels.get(lab, 0) + 1
        return {"m": r.m, "n": r.n, "n_public": r.n_public,
                "n_witness": r.n - 1 - r.n_public, "nnz": nnz, "labels": labels}


def mach_phan_tich(N=None, k=5, nonzero=False):
    """Mạch ví dụ của sách: công khai N; bí mật a, b; a, b vừa k bit; a·b = N.

    nonzero=True thêm hai cổng a - 1 ≠ 0, b - 1 ≠ 0 (cách viết tường minh của
    a, b >= 2); với N = 391 > 2^k - 1 hai cổng này thừa (xem ch05/).
    Trả về (circuit, các biến (N, a, b)).
    """
    c = Circuit()
    vN = c.public("N")
    a = c.private("a")
    b = c.private("b")
    c.to_bits(a, k, prefix="a")
    c.to_bits(b, k, prefix="b")
    if nonzero:
        c.assert_nonzero(a - 1, name="u_a")
        c.assert_nonzero(b - 1, name="u_b")
    c.assert_mul(a, b, vN)
    return c, (vN, a, b)


def mach_phan_tich_sua(N=None, k=5, bounds=(19, 23)):
    """Mạch đã sửa của Chương 12: mach_phan_tich + a <= bounds[0], b <= bounds[1].

    Với cận (19, 23): a·b <= 437 < 449 nên a·b ≡ N (mod 449) kéo theo a·b = N
    trên Z; nghiệm duy nhất của N = 391 là (17, 23). Trả về (circuit, (N, a, b)).
    """
    c = Circuit()
    vN = c.public("N")
    a = c.private("a")
    b = c.private("b")
    c.to_bits(a, k, prefix="a")
    c.to_bits(b, k, prefix="b")
    c.assert_le(a, bounds[0], prefix="da")
    c.assert_le(b, bounds[1], prefix="db")
    c.assert_mul(a, b, vN)
    return c, (vN, a, b)
