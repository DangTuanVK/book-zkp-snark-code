"""Ch.10: NTT radix-2 có đếm phép nhân: (n/2)·log2(n) phép nhân trong bướm so với n^2
của cách tính ngây thơ. Kiểm trên F_449 (n = 32, khớp snarklib.fft) và trên F_r của
BLS12-381 (n = 2^10, ω = 7^((r-1)/n); phần tử sinh 7 theo ark-bls12-381-0.6.0/src/fields/fr.rs:5)."""
import os, sys
GOC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, GOC); sys.path.insert(0, os.path.join(GOC, "ch10"))
from snarklib.field import GF
from snarklib.fft import ntt
from ho_duong_cong import ARK_R_BLS12_381


def ntt_dem(a, w, p, dem):
    """NTT đệ quy Cooley–Tukey: trả về [f(w^j)]; dem['nhan'] đếm phép nhân t·v trong bướm."""
    n = len(a)
    if n == 1:
        return list(a)
    chan = ntt_dem(a[0::2], w * w % p, p, dem)     # f_chẵn(w^2 ...)
    le = ntt_dem(a[1::2], w * w % p, p, dem)       # f_lẻ(w^2 ...)
    out, t = [0] * n, 1
    for k in range(n // 2):
        v = t * le[k] % p; dem["nhan"] += 1
        out[k], out[k + n // 2] = (chan[k] + v) % p, (chan[k] - v) % p
        t = t * w % p                              # cập nhật lũy thừa của w (không tính vào bướm)
    return out


def ngay_tho(a, w, p, dem):
    out = []
    for j in range(len(a)):
        x, s, xj = 1, 0, pow(w, j, p)
        for c in a:
            s = (s + c * x) % p; x = x * xj % p
            dem["nhan"] += 1
        out.append(s)
    return out


for ten, p, n, g in (("F_449", 449, 32, 3), ("F_r BLS12-381", ARK_R_BLS12_381, 1 << 10, 7)):
    w = pow(g, (p - 1) // n, p)
    assert pow(w, n // 2, p) == p - 1             # w có bậc đúng n
    a = [(17 * i * i + 23 * i + 391) % p for i in range(n)]
    d1, d2 = {"nhan": 0}, {"nhan": 0}
    nhanh, cham = ntt_dem(a, w, p, d1), ngay_tho(a, w, p, d2)
    assert nhanh == cham
    lg = n.bit_length() - 1
    print("%s, n = %d: NTT %d phép nhân bướm ((n/2)log n = %d), ngây thơ %d (n^2 = %d), gấp %.0f lần"
          % (ten, n, d1["nhan"], n // 2 * lg, d2["nhan"], n * n, d2["nhan"] / d1["nhan"]))
    assert d1["nhan"] == n // 2 * lg
    if p == 449:
        F = GF(449)
        assert [int(x) for x in ntt(F, a)] == nhanh and w == 221   # ω_32 = 221 (§0.4)
        print("  khớp snarklib.fft.ntt; ω_32 = %d" % w)
