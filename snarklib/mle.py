"""Mở rộng đa tuyến tính (MLE) của bảng 2^l giá trị.

QUY ƯỚC THỨ TỰ BIẾN (big-endian): chỉ số i của bảng ứng với điểm Boolean
(b_1, b_2, ..., b_l) với i = b_1·2^(l-1) + b_2·2^(l-2) + ... + b_l,
tức x_1 là bit CAO NHẤT. Ví dụ l = 2: T[0]=f(0,0), T[1]=f(0,1), T[2]=f(1,0), T[3]=f(1,1).

    f~(x) = Σ_{b ∈ {0,1}^l} T[b] · eq(b, x)
    eq(x, y) = Π_k ( x_k y_k + (1 - x_k)(1 - y_k) )
"""


def eq(F, x, y):
    """Đa thức eq(x, y) tại hai vector cùng độ dài."""
    if len(x) != len(y):
        raise ValueError("độ dài khác nhau")
    acc = F(1)
    for xk, yk in zip(x, y):
        xk, yk = F(xk), F(yk)
        acc = acc * (xk * yk + (1 - xk) * (1 - yk))
    return acc


def eq_table(F, r):
    """Bảng [eq(b, r) cho mọi b ∈ {0,1}^l] theo thứ tự big-endian, chi phí O(2^l)."""
    table = [F(1)]
    for rk in r:                     # mỗi biến mới là bit thấp hơn biến trước
        rk = F(rk)
        nxt = []
        for t in table:
            nxt.append(t * (1 - rk))   # b_k = 0
            nxt.append(t * rk)         # b_k = 1
        table = nxt
    return table


def num_vars(table):
    n = len(table)
    if n == 0 or n & (n - 1):
        raise ValueError("độ dài bảng phải là lũy thừa của 2")
    return n.bit_length() - 1


def evaluate(F, table, r):
    """Đánh giá f~(r) = Σ_i T[i] · eq_table(r)[i], chi phí O(2^l)."""
    if len(r) != num_vars(table):
        raise ValueError("số biến không khớp")
    return sum((F(t) * e for t, e in zip(table, eq_table(F, r))), F(0))


def fold(F, table, r1):
    """Cố định biến ĐẦU TIÊN x_1 = r1 (bit cao nhất), trả về bảng 2^(l-1) phần tử.

    T'[j] = (1 - r1)·T[j] + r1·T[j + 2^(l-1)].
    """
    half = len(table) // 2
    num_vars(table)
    r1 = F(r1)
    return [(1 - r1) * F(table[j]) + r1 * F(table[j + half]) for j in range(half)]


def evaluate_by_folding(F, table, r):
    """Đánh giá bằng cách fold lần lượt x_1, x_2, ..., x_l (cũng O(2^l))."""
    t = [F(x) for x in table]
    for rk in r:
        t = fold(F, t, rk)
    assert len(t) == 1
    return t[0]


def bits(i, l):
    """Điểm Boolean (b_1, ..., b_l) ứng với chỉ số i (b_1 là bit cao nhất)."""
    return [(i >> (l - 1 - k)) & 1 for k in range(l)]
