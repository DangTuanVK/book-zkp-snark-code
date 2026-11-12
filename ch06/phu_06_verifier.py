"""Chương 6 -- Kiểm chứng phụ cho Verifier (Mục 6.2, phần Verifier), không in trong sách.

  (1) Đa thức tuyến tính hóa r(X) của proof seed 2019: bậc, r(ζ) = 0, số hạng hằng
      r(0)... và đối chiếu r_0 của Verifier (bước 8) với hằng số trong r(X).
  (2) Prover gian "bỏ phần dư" với witness sai a19 = 18 (như Mã nguồn 6.4): hiệu
      tử số - t·Z_H chính là phần dư R(X), bậc < n. Đếm nghiệm của R trên F_449 và
      trong H: Verifier (với giá trị mở trung thực) chấp nhận đúng khi R(ζ) = 0.
  (3) Tỉ lệ lọt thực nghiệm của Prover gian đó qua 600 nhãn transcript.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib import plonk
from snarklib.plonk import F
from snarklib.dsl import mach_phan_tich

c = mach_phan_tich()[0]
circ = c.compile_plonkish()
a, b, cc, x = c.plonkish_assignment({"N": 391, "a": 17, "b": 23})
pk, vk = plonk.setup(circ, tau=123)
n, om = pk.n, pk.omega

P = plonk.Prover(pk, a, b, cc, x, random.Random(2019))
P.run()
r = P.linearization()
print("(1) bậc r(X) = %d; r(ζ) = %d; bậc W_ζ, W_ζω = %d, %d"
      % (r.degree(), int(r(P.zeta)), r.degree() - 1, P.z_p.degree() - 1))
tr = {}
plonk.verify(vk, x, P.proof(), trace=tr)
print("    r_0 (bước 8) = %d; v = %d, u = %d; [E]_1 = %d·[1]_1" % (int(tr["r0"]), int(tr["v"]),
      int(tr["u"]), int(tr["E"])))

# (2) Prover gian bỏ phần dư
a_sai = list(a); a_sai[19] = F(18)
Q = plonk.Prover(pk, a_sai, b, cc, x, random.Random(7), check=False)
Q.run()
_, R = Q.quotient_naive()
goc = [F(t) for t in range(449) if R(F(t)) == 0]
H = [om ** j for j in range(n)]
trong_H = [j for j in range(n) if R(H[j]) == 0]
print("(2) witness sai a19 = 18: bậc phần dư R = %d; số nghiệm của R trên F_449 = %d,"
      " trong đó %d nằm trong H" % (R.degree(), len(goc), len(trong_H)))
print("    các hàng j có R(ω^j) ≠ 0:", [j for j in range(n) if j not in trong_H])
print("    cận (n - 1)/r = %.4f; tỉ lệ nghiệm %d/449 = %.4f" % ((n - 1) / 449, len(goc), len(goc) / 449))

# (3) tỉ lệ lọt thực nghiệm
lot = lot_H = loi = 0
SO = 600
for k in range(SO):
    lab = b"plonk/lot-%d" % k
    try:
        p = plonk.Prover(pk, a_sai, b, cc, x, random.Random(1000 + k), label=lab, check=False)
        pf = p.run()
    except ZeroDivisionError:       # ζ = ω^0 = 1: L_0(ζ) không xác định, bỏ lượt
        loi += 1
        continue
    if plonk.verify(vk, x, pf, label=lab):
        lot += 1
        lot_H += p.zeta ** n == 1
print("(3) Prover gian bỏ phần dư, %d nhãn: lọt %d (%.1f%%), trong đó ζ ∈ H: %d; ζ = 1 (bỏ lượt): %d"
      % (SO, lot, 100 * lot / SO, lot_H, loi))
