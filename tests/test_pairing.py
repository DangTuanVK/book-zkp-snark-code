"""Kiểm pairing.py: song tuyến tính, không suy biến, đối xứng, e^449 = 1."""
import _path  # noqa: F401
from snarklib.curve import toy_generator
from snarklib.pairing import e, embedding_degree
from snarklib.field import GF2

E, P, G = toy_generator()
K = GF2(5387)


def test_embedding_degree():
    q, r = 5387, 449
    assert (q * q - 1) % r == 0 and (q - 1) % r != 0
    assert embedding_degree(q, r) == 2


def test_nondegenerate():
    g = e(G, G)
    assert g != 1 and g ** 449 == 1
    assert e(G, E.infinity()) == 1


def test_bilinear():
    g = e(G, G)
    pairs = [(1, 1), (2, 3), (17, 23), (28, 30), (100, 448), (448, 448), (5, 0)]
    for a, b in pairs:
        assert e(a * G, b * G) == g ** (a * b)
    assert e(G, 3 * G) * e(G, 4 * G) == e(G, 7 * G)


def test_symmetric():
    for a, b in [(2, 9), (31, 7), (400, 13)]:
        assert e(a * G, b * G) == e(b * G, a * G)


if __name__ == "__main__":
    _path.run_module(globals())
