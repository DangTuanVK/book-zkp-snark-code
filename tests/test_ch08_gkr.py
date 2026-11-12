"""Test Chương 8: GKR trên mạch phân tầng."""
import random
import _path  # noqa: F401
from snarklib.transcript import Transcript
from snarklib import gkr

rng = random.Random(81)


def test_example_circuit():
    circ, x = gkr.vi_du_mach()
    W = circ.evaluate(x)
    assert W[0] == [333, 75] and W[1] == [391, 391, 40, 58]
    proof = gkr.prove(circ, x, Transcript(b"g"))
    assert gkr.verify(circ, x, W[0], proof, Transcript(b"g"))
    assert gkr.proof_size(circ) == 34


def test_wrong_output_and_tamper():
    circ, x = gkr.vi_du_mach()
    acc = 0
    for k in range(4):
        proof = gkr.prove(circ, x, Transcript(b"w%d" % k))
        acc += gkr.verify(circ, x, [334, 75], proof, Transcript(b"w%d" % k))
        L = proof["layers"][1]
        z = (L["z"][0] + 1, L["z"][1])
        tampered = dict(proof, layers=[proof["layers"][0], dict(L, z=z)])
        acc += gkr.verify(circ, x, proof["outputs"], tampered, Transcript(b"w%d" % k))
    assert acc <= 1
    # đầu vào khác (Verifier tự tính W_d~) -> bác
    proof = gkr.prove(circ, x, Transcript(b"i"))
    assert not gkr.verify(circ, [17, 23, 28, 31], proof["outputs"], proof, Transcript(b"i"))


def test_random_circuits():
    for trial in range(5):
        sizes = [2, 4, 8, 8]
        layers = []
        for i in range(len(sizes) - 1):
            layers.append([(rng.choice(["add", "mul"]), rng.randrange(sizes[i + 1]),
                            rng.randrange(sizes[i + 1])) for _ in range(sizes[i])])
        circ = gkr.LayeredCircuit(layers, sizes[-1])
        x = [rng.randrange(449) for _ in range(sizes[-1])]
        proof = gkr.prove(circ, x, Transcript(b"r%d" % trial))
        assert gkr.verify(circ, x, circ.evaluate(x)[0], proof, Transcript(b"r%d" % trial))


if __name__ == "__main__":
    _path.run_module(globals())
