"""NTT/INTT radix-2 (Cooley–Tukey) trên F_p có 2-adicity đủ lớn.

Quy ước: ntt(F, a) trả về [f(1), f(ω), ..., f(ω^(n-1))] với f(X) = Σ a_k X^k,
ω = F.root_of_unity(n), n = len(a) là lũy thừa của 2.
Trên F_449: 449 - 1 = 2^6 · 7 nên n tối đa là 64.
"""


def _check_size(F, n):
    if n == 0 or n & (n - 1):
        raise ValueError("n phải là lũy thừa của 2")
    if (F.p - 1) % n:
        raise ValueError("F_%d không có căn bậc %d của đơn vị" % (F.p, n))


def _bit_reverse(a):
    n = len(a)
    bits = n.bit_length() - 1
    out = list(a)
    for i in range(n):
        j = int(format(i, "0%db" % bits)[::-1], 2) if bits else 0
        if i < j:
            out[i], out[j] = out[j], out[i]
    return out


def ntt_with_root(F, a, w):
    """NTT lặp (đảo bit rồi các tầng bướm) với căn đơn vị w cho trước."""
    n = len(a)
    a = _bit_reverse([F(x) for x in a])
    length = 2
    while length <= n:
        wl = w ** (n // length)          # căn bậc 'length' của đơn vị
        for start in range(0, n, length):
            t = F(1)
            for k in range(length // 2):
                u = a[start + k]
                v = a[start + k + length // 2] * t
                a[start + k] = u + v                 # bướm: u + t·v
                a[start + k + length // 2] = u - v   #       u - t·v
                t = t * wl
        length *= 2
    return a


def ntt(F, a):
    """Biến đổi thuận: hệ số -> giá trị trên H = <ω>."""
    _check_size(F, len(a))
    return ntt_with_root(F, a, F.root_of_unity(len(a)))


def intt(F, vals):
    """Biến đổi ngược: giá trị trên H -> hệ số (dùng ω^-1 và nhân n^-1)."""
    n = len(vals)
    _check_size(F, n)
    w_inv = F.root_of_unity(n).inverse()
    a = ntt_with_root(F, vals, w_inv)
    n_inv = F(n).inverse()
    return [x * n_inv for x in a]


def coset_ntt(F, a, k):
    """Đánh giá f trên coset k·H: trả về [f(k), f(kω), ..., f(kω^(n-1))].

    Thực hiện bằng cách nhân hệ số a_j với k^j rồi NTT thông thường.
    """
    k = F(k)
    scaled, kp = [], F(1)
    for x in a:
        scaled.append(F(x) * kp)
        kp = kp * k
    return ntt(F, scaled)


def coset_intt(F, vals, k):
    """Ngược của coset_ntt: giá trị trên k·H -> hệ số."""
    k_inv = F(k).inverse()
    a = intt(F, vals)
    out, kp = [], F(1)
    for x in a:
        out.append(x * kp)
        kp = kp * k_inv
    return out


def naive_eval(F, a, points):
    """Đánh giá trực tiếp (Horner) tại từng điểm — dùng để đối chiếu NTT."""
    out = []
    for x in points:
        acc = F(0)
        for c in reversed(a):
            acc = acc * x + c
        out.append(acc)
    return out
