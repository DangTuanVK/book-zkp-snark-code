"""Ch.10 (lst:ch10-tham-so): suy p, r của BN254 và BLS12-381 từ tham số u,
kiểm nguyên tố, số bit, bậc nhúng 12, 2-adicity của r - 1 và phân tích
(p^12 - 1)/r = (p^6 - 1)(p^2 + 1) · (p^4 - p^2 + 1)/r.
u, p, r đối chiếu với arkworks 0.6.0 (file:dòng trong ch10/ho_duong_cong.py)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ho_duong_cong import (bn, bls12, miller_rabin, bac_nhung, two_adicity,
                           naf_gia_tri, U_BN254, U_BLS12_381, ARK_P_BN254, ARK_R_BN254,
                           ARK_P_BLS12_381, ARK_R_BLS12_381, ARK_H1_BLS12_381,
                           ARK_ATE_LOOP_BN254)

CAC_DUONG_CONG = [("BN254", bn, U_BN254, ARK_P_BN254, ARK_R_BN254),
                  ("BLS12-381", bls12, U_BLS12_381, ARK_P_BLS12_381, ARK_R_BLS12_381)]

for ten, ho, u, p_ark, r_ark in CAC_DUONG_CONG:
    p, r, t = ho(u)                       # suy từ đa thức tham số
    n_E = p + 1 - t                       # #E(F_p) = p + 1 - t
    h = n_E // r                          # đồng hệ số (cofactor) của G1
    print("=== %s: u = %s%s (%d bit, trọng số Hamming %d)"
          % (ten, "-" if u < 0 else "", hex(abs(u)), abs(u).bit_length(), bin(abs(u)).count("1")))
    print("p =", p)
    print("  = %s (%d bit), p mod 4 = %d, p mod 3 = %d" % (hex(p), p.bit_length(), p % 4, p % 3))
    print("r =", r)
    print("  = %s (%d bit)" % (hex(r), r.bit_length()))
    print("t = %d, #E = h·r với h = %d" % (t, h))
    print("2-adicity(r - 1) = %d, 2-adicity(p - 1) = %d" % (two_adicity(r - 1), two_adicity(p - 1)))
    k = bac_nhung(p, r)
    print("bậc nhúng k = %d; r không chia hết p^j - 1 với mọi j < 12:" % k,
          all(pow(p, j, r) != 1 for j in range(1, 12)))
    # Phân tích số mũ cuối: phần dễ (p^6-1)(p^2+1), phần khó Φ12(p)/r
    phi12 = p**4 - p**2 + 1
    de, kho = (p**6 - 1) * (p**2 + 1), phi12 // r
    print("(p^12-1)/r = (p^6-1)(p^2+1)·Φ12(p)/r:", de * kho == (p**12 - 1) // r,
          "| r | Φ12(p):", phi12 % r == 0)
    print("số bit: p^12 = %d, (p^12-1)/r = %d, phần dễ = %d, phần khó = %d"
          % ((p**12).bit_length(), ((p**12 - 1) // r).bit_length(), de.bit_length(), kho.bit_length()))
    mr_p, mr_r = miller_rabin(p), miller_rabin(r)
    print("Miller–Rabin (40 cơ sở): p nguyên tố = %s, r nguyên tố = %s" % (mr_p, mr_r))
    print("khớp arkworks: p = %s, r = %s" % (p == p_ark, r == r_ark))
    assert mr_p and mr_r and k == 12 and n_E % r == 0
    assert p == p_ark and r == r_ark and p % 4 == 3 and p % 3 == 1 and de * kho == (p**12 - 1) // r
    if ten == "BN254":
        assert h == 1 and naf_gia_tri(ARK_ATE_LOOP_BN254) == 6 * u + 2
        print("NAF ATE_LOOP_COUNT của arkworks = 6u + 2: True (%d chữ số khác 0)"
              % sum(1 for d in ARK_ATE_LOOP_BN254 if d))
    else:
        assert h == (u - 1) ** 2 // 3 == ARK_H1_BLS12_381
        print("h = (u-1)^2/3 khớp COFACTOR của arkworks: True (%d bit)" % h.bit_length())
    print()
