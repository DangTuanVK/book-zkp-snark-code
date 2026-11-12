"""Chương 4 — Chi tiết phép kiểm của Ligero trên ví dụ N = 391 (Ví dụ ex:ch4-ligero).

Chạy lại đúng Mã nguồn 4.2 (cùng nhãn transcript "ch4-ligero", t = 4) và in thêm
các đại lượng trung gian mà Verifier tính: hệ số gamma của phép kiểm lân cận,
hai từ mã Enc(u_prox), Enc(u_eval) và hai vế của các đẳng thức kiểm tại từng
cột được mở. Phần cuối gập FRI trên hệ số của p(X) = 3 + X + 4X^2 + X^3 để minh
họa quan sát của Basefold (Mục 4.2c).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.field import GF
from snarklib.transcript import Transcript
from snarklib.mle import eq_table
from snarklib import ligero

F = GF(449)
N, a, b = 391, 17, 23
bits = lambda v: [(v >> i) & 1 for i in range(5)]
z = [1, N, a, b] + bits(a) + bits(b) + [0, 0]
root, st = ligero.commit(z)
r = [F(2), F(3), F(5), F(7)]
proof = ligero.open(st, r, Transcript(b"ch4-ligero"), t=4)

# Verifier dựng lại gamma từ transcript (như ligero._prelude)
tr = Transcript(b"ch4-ligero")
tr.absorb(b"root", root); tr.absorb(b"r", r); tr.absorb(b"v", F(proof["v"]))
gamma = [tr.challenge(b"gamma%d" % i) for i in range(st["m"])]
L = eq_table(F, r[:2])
print("omega_16 =", F.root_of_unity(16), " miền RS = <omega_16>, 16 phần tử")
print("gamma =", gamma, " alpha = eq(., (2,3)) =", L)
e_prox = ligero.rs_encode(proof["u_prox"])
e_eval = ligero.rs_encode(proof["u_eval"])
print("Enc(u_prox) =", e_prox)
print("Enc(u_eval) =", e_eval)
for j, col, _ in proof["cols"]:
    s_prox = sum((g * x for g, x in zip(gamma, col)), F(0))
    s_eval = sum((l * x for l, x in zip(L, col)), F(0))
    print("  cột %2d: <gamma,cột> = %3d = Enc(u_prox)[%d] = %3d;  <alpha,cột> = %3d = Enc(u_eval)[%d] = %3d"
          % (j, s_prox, j, e_prox[j], s_eval, j, e_eval[j]))
    assert s_prox == e_prox[j] and s_eval == e_eval[j]

# Hệ thống tính: 4 phần tử đầu của Enc(.) KHÔNG phải thông điệp (mã RS dạng hệ số)
print("Enc(M[1]) =", ligero.rs_encode(st["M"][1]), " (hàng (1,0,0,0) cho từ mã hằng)")

# Gập kiểu FRI trên hệ số = cố định từng biến của đa thức đa tuyến tính cùng hệ số
coef = [F(c) for c in (3, 1, 4, 1)]
fold = lambda c, al: [c[2 * i] + al * c[2 * i + 1] for i in range(len(c) // 2)]
for al1, al2 in ((F(5), F(25)), (F(2), F(3))):
    c1 = fold(coef, al1); c2 = fold(c1, al2)
    Q = coef[0] + coef[1] * al1 + coef[2] * al2 + coef[3] * al1 * al2
    print("gập (alpha1, alpha2) = (%d, %d): %s -> %s -> %s ;  Q(alpha1, alpha2) = %d"
          % (al1, al2, coef, c1, c2, Q))
    assert c2[0] == Q
