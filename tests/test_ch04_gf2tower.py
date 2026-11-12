"""Kiểm tiên đề trường cho tháp nhị phân T_1 = F_4, T_2 = F_16, T_3 = F_256 (vét cạn)."""
import _path  # noqa: F401
from snarklib.gf2tower import mul, inv, add, gen, power, fmt, T


def test_irreducible():
    """X_ι^2 + X_{ι-1} X_ι + 1 không có nghiệm trong T_ι (ι = 0, 1, 2)."""
    for lvl in range(3):
        c = gen(lvl)
        for x in range(1 << (1 << lvl)):
            assert mul(x, x, lvl) ^ mul(c, x, lvl) ^ 1 != 0


def test_f16_full_axioms():
    S = range(16)
    for a in S:
        assert mul(a, 1, 2) == a and mul(a, 0, 2) == 0
        for b in S:
            assert mul(a, b, 2) == mul(b, a, 2)
            for c in S:
                assert mul(a, mul(b, c, 2), 2) == mul(mul(a, b, 2), c, 2)
                assert mul(a, add(b, c), 2) == add(mul(a, b, 2), mul(a, c, 2))
        if a:
            assert mul(a, inv(a, 2), 2) == 1


def test_f256_axioms():
    lvl = 3
    tab = [[mul(a, b, lvl) for b in range(256)] for a in range(256)]
    for a in range(256):
        assert tab[a][1] == a
        for b in range(256):
            assert tab[a][b] == tab[b][a]
            # phân phối (vét cạn): b -> a·b là F_2-tuyến tính <=> a·b = XOR a·(bit của b)
            acc = 0
            for i in range(8):
                if b >> i & 1:
                    acc ^= tab[a][1 << i]
            assert tab[a][b] == acc
        if a:
            assert tab[a][inv(a, lvl)] == 1
    # kết hợp: tìm phần tử sinh g, kiểm g^i·g^j = g^{i+j} với mọi i, j
    g = next(x for x in range(2, 256) if len({power(x, e, lvl) for e in range(255)}) == 255)
    pw = [power(g, e, lvl) for e in range(255)]
    assert sorted(pw) == list(range(1, 256))
    for i in range(255):
        for j in range(255):
            assert tab[pw[i]][pw[j]] == pw[(i + j) % 255]


def test_subfields():
    """Tháp lồng nhau: phép nhân của T_2 trùng với phép nhân trong T_3."""
    for a in range(16):
        for b in range(16):
            assert mul(a, b, 3) == mul(a, b, 2)
    for a in range(4):
        for b in range(4):
            assert mul(a, b, 2) == mul(a, b, 1)
    # Frobenius: x^16 = x trên F_16, x^256 = x trên F_256
    assert all(power(x, 16, 2) == x for x in range(16))
    assert all(power(x, 256, 3) == x for x in range(256))


def test_hand_example():
    F16 = T(2)
    a, b = F16(0b1011), F16(0b1010)          # 1 + X0 + X0·X1,  X0 + X0·X1
    assert a * b == F16(0b1110)              # = X0 + X1 + X0·X1
    assert fmt(0b1011, 2) == "1 + X0 + X0·X1" and fmt(0b1110, 2) == "X0 + X1 + X0·X1"
    assert a * a.inverse() == F16(1) and (a / b) * b == a


if __name__ == "__main__":
    _path.run_module(globals())
