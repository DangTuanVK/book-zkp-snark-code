"""Cây Merkle nhị phân dùng SHA-256.

- Lá:      H(0x00 ‖ dữ liệu)       (tiền tố tách lá với nút trong)
- Nút trong: H(0x01 ‖ trái ‖ phải)
Số lá phải là lũy thừa của 2 (đủ cho các ví dụ trong sách).
Lá được đánh số 0..n-1 từ trái sang phải.
"""

import hashlib


def hash_leaf(data):
    if isinstance(data, str):
        data = data.encode("utf-8")
    elif isinstance(data, int):
        data = data.to_bytes(max(1, (data.bit_length() + 7) // 8), "big")
    return hashlib.sha256(b"\x00" + bytes(data)).digest()


def hash_node(left, right):
    return hashlib.sha256(b"\x01" + left + right).digest()


class MerkleTree:
    """Cây Merkle: levels[0] là các lá đã băm, levels[-1] = [gốc]."""

    def __init__(self, leaves):
        n = len(leaves)
        if n == 0 or n & (n - 1):
            raise ValueError("số lá phải là lũy thừa của 2")
        self.n = n
        self.levels = [[hash_leaf(x) for x in leaves]]
        while len(self.levels[-1]) > 1:
            cur = self.levels[-1]
            self.levels.append([hash_node(cur[2 * k], cur[2 * k + 1])
                                for k in range(len(cur) // 2)])

    @property
    def root(self):
        return self.levels[-1][0]

    def open(self, index):
        """Đường dẫn xác thực: danh sách băm anh em từ lá lên gốc."""
        if not 0 <= index < self.n:
            raise IndexError("chỉ số lá ngoài phạm vi")
        path = []
        for level in self.levels[:-1]:
            path.append(level[index ^ 1])
            index >>= 1
        return path


def commit(leaves):
    """Cam kết: trả về (gốc, cây)."""
    t = MerkleTree(leaves)
    return t.root, t


def verify(root, index, leaf, path):
    """Kiểm đường dẫn: tính lại gốc từ lá và các băm anh em."""
    h = hash_leaf(leaf)
    for sib in path:
        h = hash_node(h, sib) if index % 2 == 0 else hash_node(sib, h)
        index >>= 1
    return h == root
