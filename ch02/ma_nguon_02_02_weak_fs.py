"""Fiat--Shamir mạnh và yếu trên một PIOP đồ chơi (Chương 2).
Quan hệ R6 = {(N; a) : a(a + 6) = N trong F_449}; N = 391 có witness a = 17.
Prover cam kết A (bậc <= 2, A(1) = a) và Q = (A(A + 6) - N)/(X - 1) bằng KZG,
mở tại zeta; Verifier kiểm A(zeta)(A(zeta) + 6) - N = Q(zeta)(zeta - 1)."""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.field import GF
from snarklib.poly import Poly
from snarklib.transcript import Transcript
from snarklib import kzg

r = 449
F = GF(r)
srs = kzg.setup(3, tau=10)
X1, C6 = Poly(F, [-1, 1]), Poly(F, [6])                  # X - 1 và hằng 6


def zeta(N, cA, cQ, manh):
    t = Transcript(b"ch02/piop-r6")
    t.absorb(b"vk", srs.g1)
    if manh:                                            # FS mạnh: băm cả N
        t.absorb(b"N", N)
    t.absorb(b"A", cA); t.absorb(b"Q", cQ)
    return t.challenge(b"zeta")


def prove(N, A, Q, manh):
    cA, cQ = kzg.commit(srs, A), kzg.commit(srs, Q)
    z = zeta(N, cA, cQ, manh)
    return (cA, cQ) + kzg.open(srs, A, z) + kzg.open(srs, Q, z)


def verify(N, pi, manh):
    cA, cQ, yA, WA, yQ, WQ = pi
    z = zeta(N, cA, cQ, manh)
    return (kzg.verify(srs, cA, z, yA, WA) and kzg.verify(srs, cQ, z, yQ, WQ)
            and yA * (yA + 6) - N == yQ * (z - 1))


def co_witness(N):
    return any((a * (a + 6) - N) % r == 0 for a in range(r))


def gia_mao_yeu(A, Q, dich):
    """FS yếu: thử Qd = Q + d0 + d1 X, chọn N* SAU khi biết zeta."""
    cA = kzg.commit(srs, A)
    for delta in range(1, 20 * r):
        Qd = Q + Poly(F, [delta % r, delta // r])
        z = zeta(None, cA, kzg.commit(srs, Qd), False)
        Ns = int(A(z) * (A(z) + 6) - Qd(z) * (z - 1))   # = 391 - (d0 + d1 z)(z - 1)
        if (dich is None and not co_witness(Ns)) or Ns == dich:
            return Ns, Qd, delta


def gia_mao_manh(N, rng):
    """FS mạnh: D = A(A+6) - N - Q(X-1) có 4 nghiệm t_i, thử tới khi zeta in t."""
    thu = 0
    while True:
        t = rng.sample(range(2, r), 4); al = F(rng.randrange(1, r))
        D = Poly(F, [al * al])
        for ti in t:
            D = D * Poly(F, [-F(ti), 1])
        dis = F(36) + 4 * (F(N) + D(F(1)))              # a^2 + 6a = N + D(1)
        if not dis.is_qr():
            continue
        b = F(rng.randrange(r))
        A = Poly(F, [(dis.sqrt() - 6) / 2 - al - b, b, al])
        Q, du = (A * (A + C6) - Poly(F, [N]) - D).divmod(X1)
        assert du.is_zero() and Q.degree() <= 3
        thu += 1                                        # một lần gọi hàm băm
        if int(zeta(N, kzg.commit(srs, A), kzg.commit(srs, Q), True)) in t:
            return A, Q, thu


if __name__ == "__main__":
    rng = random.Random(2026)
    A = Poly(F, [17]) + Poly(F, [rng.randrange(r)]) * X1  # A(1) = 17
    Q = (A * (A + C6) - Poly(F, [391])) // X1
    print("trung thực N = 391: verify mạnh =", verify(391, prove(391, A, Q, True), True))
    Ns, Qd, d = gia_mao_yeu(A, Q, None)
    pi = prove(Ns, A, Qd, False)
    print("FS yếu: delta = %d, N* = %d, có witness: %s, verify yếu = %s, mạnh = %s"
          % (d, Ns, co_witness(Ns), verify(Ns, pi, False), verify(Ns, pi, True)))
    assert verify(Ns, pi, False) and not verify(Ns, pi, True)
    Ns, Qd, d = gia_mao_yeu(A, Q, 43)
    print("FS yếu, nhắm N* = 43: %d lần thử, verify yếu = %s"
          % (d, verify(43, prove(43, A, Qd, False), False)))
    A3, Q3, thu = gia_mao_manh(43, rng)
    pi = prove(43, A3, Q3, True)
    print("FS mạnh, N = 43: %d lần thử, verify mạnh = %s, có witness: %s"
          % (thu, verify(43, pi, True), co_witness(43)))
    assert verify(43, pi, True) and not co_witness(43)
