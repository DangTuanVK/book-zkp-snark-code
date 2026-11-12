"""Chương 9 — Chi phí số học không gốc (Mục 9.2a, eq:ch9-nonnative).

Mô phỏng một phép nhân trong F_q bên trong mạch R1CS trên F_r: viết a, b, k, c
thành l limb w bit, kiểm a·b = q·k + c trên số nguyên bằng (i) tích a(X)·b(X)
tại 2l-1 điểm, (ii) phương trình từng cột có carry. Đếm ràng buộc bằng
snarklib.dsl và so với công thức; phép nhân gốc chỉ tốn 1 ràng buộc.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib import folding as fo
from snarklib.field import GF

rng = random.Random(9)
F449 = GF(449)
Fp = fo.big_field(fo.PALLAS_P)
cases = [("F_5387 trong F_449", 5387, F449, w) for w in (1, 2, 3)] + \
        [("F_q(Vesta) trong F_p(Pallas)", fo.PALLAS_Q, Fp, w) for w in (32, 64, 85)]

print("trường mô phỏng          |  w |  l | wc | hợp lệ | ràng buộc | công thức | a,b đã kiểm khoảng")
for ten, q, F, w in cases:
    prm = fo.nonnative_params(q, F.p, w)
    if not prm["ok"]:                     # tràn trường: hệ số tích hoặc cột ≥ r
        print("%-24s | %2d | %2d | %2d |  không  | (cột tới %d ≥ %d)" % (ten, w, prm["l"], prm["wc"],
              prm["col_bound"], F.p))
        continue
    c, prm, assign = fo.nonnative_mul_circuit(q, F, w)
    m = c.stats()["m"]
    f_all = fo.nonnative_formula(prm["l"], w, prm["wc"])
    f_kc = fo.nonnative_formula(prm["l"], w, prm["wc"], range_ab=False)
    print("%-24s | %2d | %2d | %2d |   có   | %9d | %9d | %d" % (ten, w, prm["l"], prm["wc"], m, f_all, f_kc))
    assert m == f_all
    c2, _, _ = fo.nonnative_mul_circuit(q, F, w, range_ab=False)
    assert c2.stats()["m"] == f_kc
    # witness đúng thỏa; c sai lệch 1 thì không thỏa
    for _ in range(3):
        a, b = rng.randrange(q), rng.randrange(q)
        assert c.is_satisfied_env(c.solve(assign(a, b)))
        assert not c.is_satisfied_env(c.solve(assign(a, b, cc=(a * b) % q + 1, kk=a * b // q)))

# Ví dụ số: 1234·4321 mod 5387 trong F_449 với w = 2
a, b, q = 1234, 4321, 5387
print("Ví dụ: %d·%d = %d = 5387·%d + %d" % (a, b, a * b, a * b // q, a * b % q))
print("limb (w=2, thấp trước): a =", fo._limbs(a, 2, 7), " q =", fo._limbs(q, 2, 7))
print("Công thức: 4l(w+1) + 2(2l-1) + (2l-2)(wc+1); phép nhân gốc: 1 ràng buộc")
