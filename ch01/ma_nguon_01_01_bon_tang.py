"""Một vòng qua bốn tầng của SNARK bằng snarklib (Chương 1), mệnh đề N = 391.
Không phải SNARK hoàn chỉnh: chỉ cho thấy dữ liệu đi qua các tầng ra sao."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib import kzg
from snarklib.dsl import mach_phan_tich
from snarklib.fft import intt
from snarklib.poly import Poly, domain
from snarklib.transcript import Transcript

F = kzg.F                                    # F_449
nhap = {"N": 391, "a": 17, "b": 23}

# Tầng 1 -- arithmetization: mạch -> R1CS và Plonkish
c, _ = mach_phan_tich()
r1 = c.compile_r1cs()
x, w = c.r1cs_assignment(nhap)
print("Tầng 1: R1CS %d ràng buộc, %d biến; thỏa: %s"
      % (r1.m, r1.n, r1.is_satisfied(x, w)))
pk = c.compile_plonkish(pad_pow2=True)
a, b, o, xs = c.plonkish_assignment(nhap)
print("        Plonkish %d hàng (dùng %d); thỏa: %s"
      % (pk.n, pk.n_used, pk.is_satisfied(a, b, o, xs)))

# Tầng 2 -- PIOP: cột a thành đa thức f_a trên miền H = <ω_32>
n = len(a)
fa = Poly(F, intt(F, a))                     # f_a(ω^j) = a[j]
H = domain(F, n)
assert all(fa(H[j]) == a[j] for j in range(n))
print("Tầng 2: n = %d, ω_%d = %d, bậc f_a = %d"
      % (n, n, int(H[1]), fa.degree()))

# Tầng 3 -- PCS: cam kết KZG của f_a (SRS toy, cửa sập τ = 123)
srs = kzg.setup(n - 1, 123)
C = kzg.commit(srs, fa)
print("Tầng 3: C_a = [f_a(τ)]_1 = (%d, %d)" % (int(C.x), int(C.y)))

# Tầng 4 -- Fiat--Shamir: thách thức ζ rút từ transcript chứa N, C_a
tr = Transcript(b"ch01/bon-tang")
tr.absorb("N", 391)
tr.absorb("C_a", C)
zeta = tr.challenge("zeta")
y, W = kzg.open(srs, fa, zeta)
print("Tầng 4: ζ = %d, f_a(ζ) = %d, kiểm mở KZG: %s"
      % (int(zeta), int(y), kzg.verify(srs, C, zeta, y, W)))
