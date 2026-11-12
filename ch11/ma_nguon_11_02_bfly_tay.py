"""Ch.11: mô hình của ntt_bfly.sv trên F_449 (W = 16) và kiểm một NTT n = 16
dựng hoàn toàn từ butterfly này, đối chiếu snarklib.fft.ntt."""
import os, sys
GOC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, GOC)
from snarklib.montgomery import Montgomery
from snarklib.field import GF
from snarklib import fft

P, W = 449, 16
M = Montgomery(P, word=W)
F = GF(P)


def bfly(u, v, w_m):
    """Ba tầng pipeline: montmul (2 chu kỳ), sau đó cộng/trừ mod p (1 chu kỳ)."""
    t = M.mul(v, w_m)                      # t = REDC(v·w~) = v·ω mod p
    x = u + t - P if u + t >= P else u + t
    y = u - t if u >= t else u - t + P
    return t, x, y


om = int(F.root_of_unity(16))              # ω_16 = 349
w_m = om * M.R % P                         # twiddle lưu ở dạng Montgomery
u, v = 17, 23
t, x, y = bfly(u, v, w_m)
print("ω_16 = %d, w_m = ω·2^16 mod p = %d (0x%04x)" % (om, w_m, w_m))
print("u = %d, v = %d: t = ω·v mod p = %d, x = u + t = %d, y = u - t = %d" % (u, v, t, x, y))
assert t == om * v % P and x == (u + om * v) % P and y == (u - om * v) % P

# NTT n = 16: đảo bit, log2(n) = 4 tầng, mỗi tầng n/2 = 8 butterfly (tổng 32)
n = 16
a = [(3 * i * i + 7 * i + 1) % P for i in range(n)]   # hệ số f(X) bất kỳ
A = fft._bit_reverse(a)
so_bfly, length = 0, 2
while length <= n:
    wl = pow(om, n // length, P)
    for s in range(0, n, length):
        for k in range(length // 2):
            i, j = s + k, s + k + length // 2
            _, A[i], A[j] = bfly(A[i], A[j], pow(wl, k, P) * M.R % P)
            so_bfly += 1
    length *= 2
ref = [int(z) for z in fft.ntt(F, [F(z) for z in a])]
print("NTT n = %d bằng %d butterfly: %s" % (n, so_bfly, A))
print("khớp snarklib.fft.ntt: %s" % (A == ref))
assert A == ref and so_bfly == (n // 2) * 4
