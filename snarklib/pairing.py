"""Pairing Tate rút gọn (Type 1) trên nhóm con bậc r = 449 của E: y^2 = x^3 + x / F_5387.

    e(P, Q) = f_{r,P}( ψ(Q) ) ^ ((q^2 - 1) / r)

- f_{r,P}: hàm Miller, div(f) = r(P) - r(O).
- ψ(x, y) = (-x, i·y): ánh xạ méo (distortion map) đưa Q ∈ E(F_q) sang
  E(F_{q^2}) ngoài nhóm <P>, nhờ đó e(G, G) ≠ 1 (pairing đối xứng).
- Bậc nhúng k = 2: r | q^2 - 1 nhưng r ∤ q - 1.

Khử mẫu số (denominator elimination): đường thẳng đứng x - x_T đánh giá tại
ψ(Q) cho giá trị -x_Q - x_T ∈ F_q; mọi phần tử của F_q^* bị lũy thừa cuối
đưa về 1 (vì q - 1 chia hết (q^2 - 1)/r), nên ta bỏ qua các mẫu số đó.
Chỉ dùng để tính tay: bậc nhúng 2 ở kích thước thật không an toàn (MOV).
"""

from .field import GF2
from .curve import toy_curve, R_TOY


def distortion(Q, K):
    """ψ(x, y) = (-x, i·y), trả về cặp tọa độ trong F_{q^2}."""
    return K(-int(Q.x), 0), K(0, int(Q.y))


def _line(T, S, xq, yq, K):
    """Giá trị tại (xq, yq) của đường thẳng qua T và S (tiếp tuyến nếu T = S).

    l(x, y) = (y - y_T) - λ (x - x_T). Không gọi với T = -S (đường thẳng đứng).
    """
    if T == S:
        lam = (3 * T.x ** 2 + T.E.A) / (2 * T.y)
    else:
        lam = (S.y - T.y) / (S.x - T.x)
    return (yq - K(int(T.y))) - K(int(lam)) * (xq - K(int(T.x)))


def miller(P, Q, r=R_TOY):
    """Vòng lặp Miller: f_{r,P}(ψ(Q)) ∈ F_{q^2} (đã khử các mẫu số đường đứng)."""
    K = GF2(P.E.q)
    xq, yq = distortion(Q, K)
    f = K(1)
    T = P
    for bit in bin(r)[3:]:                      # bỏ bit cao nhất
        # bước nhân đôi: f <- f^2 · l_{T,T}(ψQ), T <- 2T
        if T.y == 0:
            T = T.E.infinity()                  # tiếp tuyến đứng: bỏ qua mẫu
        else:
            f = f * f * _line(T, T, xq, yq, K)
            T = T + T
        if bit == "1":
            # bước cộng: f <- f · l_{T,P}(ψQ), T <- T + P
            if T == -P:
                T = T.E.infinity()              # đường đứng (bước cuối): bỏ qua
            else:
                f = f * _line(T, P, xq, yq, K)
                T = T + P
    assert T.is_infinity(), "r·P phải bằng O"
    return f


def tate(P, Q, r=R_TOY):
    """Pairing Tate rút gọn e(P, Q) ∈ μ_r ⊂ F_{q^2}^*."""
    K = GF2(P.E.q)
    if P.is_infinity() or Q.is_infinity():
        return K(1)
    q = P.E.q
    assert (q * q - 1) % r == 0
    return miller(P, Q, r) ** ((q * q - 1) // r)


def e(P, Q):
    """Tên ngắn: e(P, Q) = tate(P, Q) với r = 449."""
    return tate(P, Q, R_TOY)


def embedding_degree(q, r, kmax=24):
    """Bậc nhúng: k nhỏ nhất với r | q^k - 1."""
    for k in range(1, kmax + 1):
        if (pow(q, k, r) - 1) % r == 0:
            return k
    return None
