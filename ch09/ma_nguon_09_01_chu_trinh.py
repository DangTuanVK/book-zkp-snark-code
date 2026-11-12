"""Chương 9 — Chu trình đường cong và 2-chain thu nhỏ (Mục 9.2a).

(1) Tìm chu trình nhỏ nhất E1/F_q1: y^2 = x^3 + b1, E2/F_q2: y^2 = x^3 + b2 với
#E1 = q2, #E2 = q1 (dạng j = 0 như Pasta); đếm điểm lại bằng vét cạn (x, y).
(2) Vì sao đường cong pairing của sách (q = 5387, r = 449) không nằm trong chu trình.
(3) Hai 2-chain chứa đường cong của sách: làm đường NGOÀI (đường trong trên F_449)
và làm đường TRONG (đường ngoài siêu kỳ dị có nhóm con cấp 5387).
"""
import math, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib import folding as fo
from snarklib.field import factorize

def dem_vet_can(p, a, b):                  # đếm điểm theo định nghĩa, O(p^2)
    return 1 + sum(1 for x in range(p) for y in range(p) if (y*y - x**3 - a*x - b) % p == 0)

# (1) Chu trình nhỏ nhất
cyc = fo.find_cycles(5, 400, b_max=10, limit=3)
for c in cyc:
    n1, n2 = dem_vet_can(c["q1"], 0, c["b1"]), dem_vet_can(c["q2"], 0, c["b2"])
    print("Chu trình: E1/F_%d: y^2=x^3+%d, #E1=%d | E2/F_%d: y^2=x^3+%d, #E2=%d | bậc nhúng %d, %d"
          % (c["q1"], c["b1"], n1, c["q2"], c["b2"], n2, c["k1"], c["k2"]))
    assert n1 == c["q2"] and n2 == c["q1"]
c = cyc[0]
G1 = fo.first_point(c["q1"], 0, c["b1"], 1)
assert fo.ec_mul(c["q2"], G1, 0, c["q1"]) is None   # mọi điểm khác O có cấp q2
print("Điểm", G1, "trên E1 có cấp", fo.point_order(G1, 0, c["q1"], c["q2"]))

# (2) Đường cong của sách: #E = 5388 = 12·449; cần đường trên F_449 có đúng 5387 điểm
q, r = 5387, 449
lo, hi = math.ceil(r + 1 - 2 * r ** 0.5), math.floor(r + 1 + 2 * r ** 0.5)
print("Hasse: #E'(F_449) ∈ [%d, %d] (|#E' - 450| ≤ 2√449 ≈ 42,4) -> không thể bằng %d" % (lo, hi, q))
print("Chu trình cần q ≈ r; ở đây q/r ≈ %d, và #E = 12·449 không nguyên tố" % round(q / r))
assert not (lo <= q <= hi)

# (3a) Đường của sách làm đường NGOÀI: đường trong trên F_449 (trường vô hướng của E)
inn = fo.chain_inner_over(r, kmax=6)
print("2-chain (E ngoài): E_in/F_449: y^2=x^3+%dx+%d, #E_in=%d=%s, nhóm con r'=%d, bậc nhúng %d"
      % (inn["a"], inn["b"], inn["n"], factorize(inn["n"]), inn["r"], inn["k"]))
assert dem_vet_can(r, inn["a"], inn["b"]) == inn["n"] and inn["n"] % inn["r"] == 0

# (3b) Đường của sách làm đường TRONG: E_out/F_P siêu kỳ dị, 5387 | P + 1
out = fo.chain_outer_supersingular(q)
print("2-chain (E trong): E_out/F_%d: y^2=x^3+x, P = %d·5387 - 1, #E_out=%d, điểm %s cấp %d, bậc nhúng %d"
      % (out["P"], out["h"], out["n"], out["G"], out["order"], out["k"]))
assert out["n"] == out["P"] + 1 and out["order"] == q and out["k"] == 2

# Pasta (số thật): Pallas/F_p có cấp q, Vesta/F_q có cấp p (halo2curves pasta)
print("Pasta: p =", hex(fo.PALLAS_P), "\n       q =", hex(fo.PALLAS_Q),
      "\n       q - p =", fo.PALLAS_Q - fo.PALLAS_P, "(≈ 2^%d)" % (fo.PALLAS_Q - fo.PALLAS_P).bit_length())
assert fo.probable_prime(fo.PALLAS_P) and fo.probable_prime(fo.PALLAS_Q)
print("OK: mọi phép đếm điểm khớp")
