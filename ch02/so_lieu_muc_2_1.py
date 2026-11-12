"""Chương 2, Mục 2.1 -- số liệu cho phần định nghĩa và mô hình (tệp phụ, không in).

(1) Quan hệ của mạch R_C so với quan hệ dự định R_N: witness và ngôn ngữ.
(2) Ví dụ biểu diễn đại số trong AGM với SRS tau = 10 (cùng ví dụ Chương 3).
(3) Trên nhóm toy, q-DLOG và q-SDH đều giải được bằng vét cạn.
(4) Cận nhóm tổng quát của Boneh--Shoup (Định lý 16.3) tại q = 449.
"""
import os, sys, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.field import GF
from snarklib.curve import toy_generator
from snarklib.pairing import e
from snarklib.dsl import mach_phan_tich

r = 449
F = GF(r)

print("== (1) R_N (số nguyên) và R_C (mạch trên F_449)")
c = mach_phan_tich()[0]
r1cs = c.compile_r1cs()


def witness_mach(N):
    """Mọi (a, b) in [0,31]^2 thỏa R1CS của mạch với đầu vào công khai N."""
    out = []
    for a in range(32):
        for b in range(32):
            x, w = c.r1cs_assignment({"N": N, "a": a, "b": b})
            if r1cs.is_satisfied(x, w):
                out.append((a, b))
    return out


for N in (391, 389):
    wn = [(a, b) for a in range(2, 32) for b in range(2, 32) if a * b == N]
    wc = witness_mach(N)
    print("N = %d: witness R_N = %s; witness R_C = %s" % (N, wn, wc))
    assert all((a * b - N) % r == 0 for a, b in wc)
print("389 + 449 = %d = %s" % (389 + 449, "2·419"))
assert 838 == 2 * 419

# ngôn ngữ: kiểm bằng số học, đối chiếu R1CS ở hai điểm trên
LN = sorted({a * b for a in range(2, 32) for b in range(2, 32)})
LC = sorted({(a * b) % r for a in range(32) for b in range(32)})
print("|L(R_N)| = %d (các tích a·b với 2 <= a,b <= 31, từ %d tới %d)"
      % (len(LN), LN[0], LN[-1]))
print("|L(R_C)| = %d / %d phần tử của F_449" % (len(LC), r))
thieu = [x for x in range(r) if x not in set(LC)]
print("phần tử không có witness trong R_C:", len(thieu), "ví dụ", thieu[:8])
dem = {}
for a in range(32):
    for b in range(32):
        dem[(a * b) % r] = dem.get((a * b) % r, 0) + 1
print("số witness R_C của N = 391:", dem[391], "; N = 0:", dem[0],
      "; trung bình trên L(R_C): %.2f" % (1024 / len(LC)))
lech = [x for x in range(4, r) if x in LN and
        dem.get(x, 0) != sum(1 for a in range(2, 32) for b in range(2, 32) if a * b == x)]
print("số N in [4,448] thuộc L(R_N) mà số witness R_C khác R_N:", len(lech))

print("\n== (2) Biểu diễn đại số (AGM), SRS tau = 10")
E, P, G = toy_generator()
tau = F(10)
srs = [G * int(tau ** i) for i in range(4)]
coef = [23, 17, 0, 1]                      # f(X) = X^3 + 17X + 23
C = srs[0] * 0
for ci, Si in zip(coef, srs):
    C = C + Si * ci
ftau = sum(F(ci) * tau ** i for i, ci in enumerate(coef))
print("G1 =", G, "; [tau^i]_1 =", srs)
print("f(tau) =", int(ftau), "; C =", C, "; C == f(tau)·G1:", C == G * int(ftau))
alt = [int(ftau), 0, 0, 0]
C2 = srs[0] * alt[0]
print("biểu diễn thứ hai", alt, "cho cùng C:", C2 == C)
# R(X) = f(X) - f(tau): hiệu hai biểu diễn, triệt tiêu tại tau
Rc = [(coef[0] - int(ftau)) % r, coef[1], coef[2], coef[3]]
roots = [t for t in range(r) if sum(F(ci) * F(t) ** i for i, ci in enumerate(Rc)) == 0]
print("R(X) = X^3 + 17X + %d (mod 449); nghiệm trong F_449: %s" % (Rc[0], roots))

print("\n== (3) Trên nhóm toy: q-DLOG và q-SDH giải bằng vét cạn")
T1 = srs[1]
buoc, Q = 0, E.infinity()
while Q != T1:
    Q = Q + G
    buoc += 1
print("vét cạn [tau]_1: tau =", buoc, "sau", buoc, "phép cộng điểm")
cc = 1
inv = F(1) / (F(buoc) + cc)
sol = G * int(inv)
lhs = e(sol, srs[1] + G * cc)
rhs = e(G, G)
print("lời giải q-SDH với c = %d: 1/(tau+c) = %d, [1/(tau+c)]_1 = %s" % (cc, int(inv), sol))
print("kiểm e([1/(tau+c)]_1, [tau]_1 + c·G1) == e(G1, G1):", lhs == rhs, "; e(G1,G1) =", rhs)

print("\n== (4) Cận nhóm tổng quát Pr <= ((3T+2)^2/2 + 1)/q")
def can(T, q):
    return ((3 * T + 2) ** 2 / 2 + 1) / q
for T in (5, 6, 7, 10, 14):
    print("q = 449, T = %2d: cận = %.4f" % (T, can(T, r)))
T_half = next(T for T in range(1, 1000) if can(T, r) >= 0.5)
T_one = next(T for T in range(1, 1000) if can(T, r) >= 1)
print("T nhỏ nhất để cận >= 1/2:", T_half, "; >= 1:", T_one)
print("baby-step/giant-step: khoảng 2·sqrt(449) = %.1f phép toán nhóm" % (2 * math.sqrt(r)))
rbls = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
# T sao cho (3T)^2/2 ~ r/2  =>  T ~ sqrt(r)/3
print("BLS12-381: log2 r = %.2f; cận >= 1/2 cần T >= khoảng 2^%.2f"
      % (math.log2(rbls), math.log2(math.sqrt(rbls) / 3)))
print("r - 1 =", "2^6·7", "; Cheon/Boneh--Shoup 16.1.4 với d = 64: sqrt(449/64) + sqrt(64) = %.2f"
      % (math.sqrt(r / 64) + 8))
