"""Chương 8 — HyperPlonk tối giản: zerocheck cổng + productcheck hoán vị (Mục 8.2d).

Mạch Plonkish của ví dụ N = 391 (20 hàng, snarklib.dsl) đặt trên siêu lập phương
{0,1}^5 (đệm 32 hàng) thay vì miền nhân bậc 32 của PLONK — không cần FFT.
  zerocheck:    Σ_x eq(τ, x)·gate(x) = 0,   gate = qL·a + qR·b + qO·c + qM·a·b + (qC + PI)
  productcheck: Π_x f1(x)/f2(x) = 1,  f1 = Π_col (w + β·id + γ),  f2 = Π_col (w + β·σ + γ)
Oracle là chính các bảng (chưa gắn PCS) — mức PIOP.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.field import GF
from snarklib.transcript import Transcript
from snarklib.dsl import mach_phan_tich
from snarklib import sumcheck as sc

F = GF(449)
circ, _ = mach_phan_tich()
pk = circ.compile_plonkish()
a, b, c, x = circ.plonkish_assignment({"N": 391, "a": 17, "b": 23})
T = sc.plonkish_tables(F, pk, a, b, c, x)
mu = len(T["a"]).bit_length() - 1
print("Plonkish: %d hàng, %d ràng buộc sao chép -> siêu lập phương µ = %d (%d hàng)"
      % (pk.n, len(pk.copies), mu, len(T["a"])))

proof, v = sc.hyperplonk_prove(F, T, Transcript(b"ch8-hp"))
ok = sc.hyperplonk_verify(F, T, proof, v, Transcript(b"ch8-hp"))
print("Zerocheck cổng: %d vòng, bậc %d mỗi vòng (gate bậc 3 nhân eq)" % (mu, len(proof["gate"]["polys"][0]) - 1))
print("Productcheck: ṽ có %d biến; tích s = ṽ(1,..,1,0) = %d; zerocheck gộp %d vòng, bậc %d"
      % (mu + 1, proof["perm"]["s"], len(proof["perm"]["polys"]), len(proof["perm"]["polys"][0]) - 1))
n_field = sum(len(s) for s in proof["gate"]["polys"]) + sum(len(s) for s in proof["perm"]["polys"]) + 1
print("Phần tử trường do sum-check gửi: %d; nonce chọn (β, γ): %d; Verifier chấp nhận: %s"
      % (n_field, proof["nonce"], ok))

# (1) Witness sai (17, 24): cổng nhân a·b = N hỏng -> zerocheck bác
a2, b2, c2, _ = circ.plonkish_assignment({"N": 391, "a": 17, "b": 24})
T2 = sc.plonkish_tables(F, pk, a2, b2, c2, x)
bad_gate = [i for i, e in enumerate(pk.gate_residuals(a2, b2, c2, x)) if e != 0]
# (2) Giữ mọi cổng đúng nhưng phá một ràng buộc sao chép: hàng nhân đầu tiên, a += 1, c tính lại
j = next(i for i in range(pk.n) if pk.sel["qM"][i] != 0 and pk.sel["qL"][i] == 0 and pk.sel["qC"][i] == 0)
a3, c3 = list(a), list(c)
a3[j] = F(a[j]) + 1
c3[j] = -(pk.sel["qM"][j] * a3[j] * b[j]) / pk.sel["qO"][j]
T3 = sc.plonkish_tables(F, pk, a3, b, c3, x)
print("Witness (17,24): cổng hỏng ở hàng", bad_gate, "; bản phá sao chép ở hàng", j,
      ": cổng còn đúng =", all(e == 0 for e in pk.gate_residuals(a3, b, c3, x)),
      ", vi phạm sao chép =", len(pk.copy_violations(a3, b, c3)))
res = []
for TT in (T2, T3):                    # 4 nhãn transcript, cho phép lọt tối đa 1 lần (~1/449)
    res.append(sum(sc.hyperplonk_verify(F, TT, *sc.hyperplonk_prove(F, TT, Transcript(b"h%d" % k)),
                                        Transcript(b"h%d" % k)) for k in range(4)))
print("Lọt qua Verifier: witness sai %d/4, phá sao chép %d/4" % tuple(res))
assert ok and pk.n == 20 and mu == 5 and res[0] <= 1 and res[1] <= 1
