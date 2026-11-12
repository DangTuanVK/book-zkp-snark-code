"""Chương 7 -- Điều kiện đặc số của LogUp trên F_449.

Đồng nhất thức Σ 1/(X + f_i) = Σ m_j/(X + t_j) chỉ "thấy" bội mod p. Một giá
trị ngoài bảng lặp đúng p = 449 lần có hệ số cực điểm 449 ≡ 0 nên biến mất:
LogUp chấp nhận tại mọi λ, trong khi kiểm tích (số mũ nguyên) vẫn bác.
Điều kiện an toàn: số phần tử tra bảng |f| < p (chờ đối chiếu Haböck 2022/1530).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib import logup
from snarklib.logup import F

# (1) phản ví dụ: f = (17, 23) ‖ 40 lặp 449 lần, bảng {0..31}
kq = logup.characteristic_counterexample(honest=(17, 23), bad=40)
print("|f| = %d (2 giá trị hợp lệ + %d lần giá trị %d ∉ bảng)"
      % (kq["f_len"], kq["outside_count"], kq["bad"]))
print("hệ số cực điểm khác 0 (vế trái - vế phải):", kq["poles"])
print("LogUp đúng tại %d/%d giá trị λ thử:" % (kq["lams"], kq["lams"]), kq["logup_ok"])
print("kiểm tích Π(λ + f_i) = Π(λ + t_j)^{m_j} đúng tại ít nhất một λ đó:", kq["product_ok_any"])
print("  (bậc hai vế của kiểm tích: %d và %d)" % (kq["f_len"], 2))

# (2) bớt một phần tử: 448 lần -> hệ số 448 ≡ -1 ≠ 0, LogUp bác
f448 = [F(17), F(23)] + [F(40)] * 448
t = logup.range_table(31)
m, _ = logup.field_multiplicities(f448, t)
print("448 lần: hệ số cực điểm", {k: int(v) for k, v in logup.pole_coefficients(f448, t, m).items()},
      "-> LogUp tại λ = 5:", logup.identity_at(f448, t, m, 5))

# (3) vì sao không dựng được trong PIOP trên F_449: cần |f| ≥ p nhưng |H| ≤ 2^6
v2 = 0
while (448 >> v2) % 2 == 0:
    v2 += 1
print("cần ≥ %d lookup; F_449 có 449 - 1 = 2^%d·7 nên miền NTT lớn nhất |H| = %d"
      % (logup.min_cheat_size(), v2, 2 ** v2))

assert kq["logup_ok"] and not kq["product_ok_any"] and kq["poles"] == {}
assert not logup.identity_at(f448, t, m, 5) and 2 ** v2 == 64 < 449
