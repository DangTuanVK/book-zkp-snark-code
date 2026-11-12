"""snarklib: thư viện Python thuần dùng chung cho sách "zk-SNARK chuyên sâu".

Chỉ dùng thư viện chuẩn để bạn đọc theo được từng phép tính.
Hằng số xuyên suốt: r = 449, q = 5387, E: y^2 = x^3 + x, N = 391.
"""

from .field import GF, GF2, is_prime, factorize

R = 449          # bậc trường vô hướng F_r
Q = 5387         # trường cơ sở của đường cong toy
N_EXAMPLE = 391  # ví dụ xuyên suốt: 17 · 23

__all__ = ["GF", "GF2", "is_prime", "factorize", "R", "Q", "N_EXAMPLE"]
