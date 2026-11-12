"""Tính mọi hằng số xuyên suốt của sách bằng chính snarklib và ghi ket_qua/hang_so.txt.

Chạy:  python -m snarklib.hang_so      (từ thư mục gốc kho mã)
Mỗi hằng số được tính lại từ đầu, không chép số có sẵn; các giá trị kỳ vọng
trong plan.md §0.4 chỉ dùng để đối chiếu (dòng "# KHÁC KỲ VỌNG" nếu lệch).
"""

import os

from .field import GF, GF2, is_prime, factorize
from .curve import toy_curve
from .pairing import e, embedding_degree

R, Q, N = 449, 5387, 391
NS = (8, 16, 32, 64)


def _fmt_fact(d):
    return "·".join("%d^%d" % (p, k) if k > 1 else "%d" % p for p, k in sorted(d.items()))


def coset_constants(F, n):
    """k1 = số nhỏ nhất ≥ 2 với k1^n ≠ 1; k2 = số nhỏ nhất > k1 với k2^n ≠ 1 và (k2/k1)^n ≠ 1."""
    k1 = next(k for k in range(2, F.p) if F(k) ** n != 1)
    k2 = next(k for k in range(k1 + 1, F.p)
              if F(k) ** n != 1 and (F(k) / F(k1)) ** n != 1)
    # kiểm H, k1·H, k2·H rời nhau trực tiếp
    w = F.root_of_unity(n)
    H = {int(w ** j) for j in range(n)}
    k1H = {int(F(k1) * h) for h in H}
    k2H = {int(F(k2) * h) for h in H}
    assert not (H & k1H) and not (H & k2H) and not (k1H & k2H)
    return k1, k2


def compute():
    """Trả về (danh sách dòng, danh sách cảnh báo lệch kỳ vọng)."""
    lines, issues = [], []

    def put(name, value, expected=None):
        lines.append("%s = %s" % (name, value))
        if expected is not None and str(value) != str(expected):
            lines.append("# KHÁC KỲ VỌNG: %s kỳ vọng %s" % (name, expected))
            issues.append("%s = %s (kỳ vọng %s)" % (name, value, expected))

    def note(s):
        lines.append("# " + s)

    F = GF(R)
    note("=== Trường vô hướng F_r ===")
    note("is_prime bằng chia thử; phân tích r-1 bằng chia thử")
    put("r", R)
    put("r_nguyen_to", is_prime(R), True)
    put("r_minus_1", "%d = %s" % (R - 1, _fmt_fact(factorize(R - 1))), "448 = 2^6·7")
    put("two_adicity_r", F.two_adicity(), 6)
    note("g = số nhỏ nhất có g^((r-1)/p) ≠ 1 với mọi ước nguyên tố p của r-1")
    g = F.generator()
    put("g", g, 3)
    note("omega_n = g^((r-1)/n); bac_omega_n = bậc nhân tính bằng FieldElement.order()")
    for n in NS:
        w = F.root_of_unity(n)
        put("omega_%d" % n, int(w), 349 if n == 16 else None)
        put("omega_%d_mu" % n, "g^%d" % ((R - 1) // n))
        put("bac_omega_%d" % n, w.order(), n)

    note("=== Đường cong toy E: y^2 = x^3 + x trên F_q ===")
    E = toy_curve()
    put("q", Q)
    put("q_nguyen_to", is_prime(Q), True)
    put("q_mod_4", Q % 4, 3)
    note("#E đếm bằng 1 + Σ_x (1 + Legendre(x^3 + x)) trên mọi x ∈ F_q")
    nE = E.count_points()
    put("so_diem_E", nE, 5388)
    put("so_diem_E_phan_tich", _fmt_fact(factorize(nE)), "2^2·3·449")
    put("cofactor", nE // R, 12)
    note("P: x nhỏ nhất với x^3+x là thặng dư bậc hai khác 0, y = căn nhỏ hơn; G1 = 12·P (≠ O)")
    P, G = E.deterministic_generator(R, nE // R)
    put("P", P)
    put("G1", G)
    put("G1_khac_O", not G.is_infinity(), True)
    put("r_G1_bang_O", (R * G).is_infinity(), True)
    put("bac_G1", G.order(R), R)
    note("bậc nhúng: k nhỏ nhất với r | q^k - 1")
    put("r_chia_q_minus_1", (Q - 1) % R == 0, False)
    put("r_chia_q2_minus_1", (Q * Q - 1) % R == 0, True)
    put("bac_nhung", embedding_degree(Q, R), 2)
    note("F_{q^2} = F_q[i]/(i^2+1); e = Tate rút gọn, Miller f_{r,P}(ψ(Q))^((q^2-1)/r), ψ(x,y)=(-x, i·y)")
    put("mu_luy_thua_cuoi", "(q^2-1)/r = %d" % ((Q * Q - 1) // R))
    eGG = e(G, G)
    put("e_G1_G1", eGG)
    put("e_G1_G1_khac_1", eGG != 1, True)
    put("e_G1_G1_mu_r_bang_1", eGG ** R == 1, True)
    K = GF2(Q)
    assert isinstance(eGG, K)

    note("=== Hằng số coset PLONK trên F_r ===")
    note("k1 = số nhỏ nhất ≥ 2 với k1^n ≠ 1; k2 = số nhỏ nhất > k1 với k2^n ≠ 1 và (k2/k1)^n ≠ 1")
    for n in NS:
        k1, k2 = coset_constants(F, n)
        put("k1_n%d" % n, k1)
        put("k2_n%d" % n, k2)

    note("=== Ví dụ xuyên suốt N = 391 ===")
    put("N", "%d = %s" % (N, _fmt_fact(factorize(N))), "391 = 17·23")
    note("duyệt mọi (a, b) ∈ [2,31]^2 với a·b ≡ N (mod r)")
    sols = [(a, b) for a in range(2, 32) for b in range(2, 32) if (a * b - N) % R == 0]
    put("nghiem_mod_r_trong_2_31", sols, [(17, 23), (23, 17), (28, 30), (30, 28)])
    put("nghiem_gia_tich_nguyen", "28·30 = %d ≡ %d (mod %d)" % (28 * 30, 28 * 30 % R, R))
    note("bản sửa: 2 ≤ a ≤ 19, 2 ≤ b ≤ 23")
    sols2 = [(a, b) for a in range(2, 20) for b in range(2, 24) if (a * b - N) % R == 0]
    put("nghiem_ban_sua", sols2, [(17, 23)])
    put("max_ab_ban_sua", 19 * 23, 437)
    put("max_ab_nho_hon_r", 19 * 23 < R, True)
    return lines, issues


def write(path=None):
    if path is None:
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        path = os.path.join(root, "ket_qua", "hang_so.txt")
    lines, issues = compute()
    header = ["# Hằng số xuyên suốt sách zk-SNARK chuyên sâu — sinh bởi snarklib/hang_so.py",
              "# Không sửa tay; chạy lại: python -m snarklib.hang_so"]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(header + lines) + "\n")
    return path, lines, issues


if __name__ == "__main__":
    p, lines, issues = write()
    print("\n".join(lines))
    print("\nĐã ghi", p)
    if issues:
        print("LỆCH KỲ VỌNG:", *issues, sep="\n  ")
