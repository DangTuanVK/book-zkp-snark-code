"""Kiểm merkle.py: commit, open, verify, phát hiện giả mạo."""
import _path  # noqa: F401
from snarklib.merkle import commit, verify, MerkleTree


def test_merkle():
    leaves = [("la-%d" % k).encode() for k in range(16)]
    root, tree = commit(leaves)
    for k in range(16):
        path = tree.open(k)
        assert len(path) == 4
        assert verify(root, k, leaves[k], path)
        assert not verify(root, k, b"sai", path)
        assert not verify(root, (k + 1) % 16, leaves[k], path)
    try:
        MerkleTree([b"a", b"b", b"c"])
        assert False
    except ValueError:
        pass


if __name__ == "__main__":
    _path.run_module(globals())
