"""Transcript Fiat–Shamir dựa trên SHA-256 có tách miền (domain separation).

    t = Transcript(b"vi-du")
    t.absorb(b"C", diem_cam_ket)
    alpha = t.challenge(b"alpha")      # phần tử của F_r, r = 449 mặc định

Trạng thái là một chuỗi băm 32 byte. Mỗi lần absorb/challenge, trạng thái
mới = SHA256(trạng thái cũ ‖ thẻ thao tác ‖ độ dài‖nhãn ‖ độ dài‖dữ liệu),
nên thứ tự, nhãn và ranh giới dữ liệu đều ảnh hưởng tới thách thức.

Thách thức: lấy 64 byte = SHA256(s‖"C0") ‖ SHA256(s‖"C1"), đọc thành số
nguyên 512 bit rồi lấy mod r. Độ lệch thống kê so với phân bố đều trên F_r
nhỏ hơn r / 2^512 (với r = 449 là < 2^-503), không đáng kể.
"""

import hashlib

from .field import GF, FieldElement, Fq2Element

DOMAIN = b"snarklib/transcript/v1"


def _u32(n):
    return n.to_bytes(4, "big")


def encode(data):
    """Mã hóa có thẻ kiểu để hai đầu vào khác kiểu không trùng chuỗi byte."""
    if isinstance(data, (bytes, bytearray)):
        return b"B" + bytes(data)
    if isinstance(data, str):
        return b"S" + data.encode("utf-8")
    if isinstance(data, bool):
        raise TypeError("không nhận bool")
    if isinstance(data, int):
        if data < 0:
            raise ValueError("chỉ nhận số nguyên không âm")
        return b"I" + data.to_bytes(max(1, (data.bit_length() + 7) // 8), "big")
    if isinstance(data, FieldElement):
        width = (data.p.bit_length() + 7) // 8
        return b"F" + _u32(data.p) + data.v.to_bytes(width, "big")
    if isinstance(data, Fq2Element):
        width = (data.q.bit_length() + 7) // 8
        return b"K" + _u32(data.q) + data.a.to_bytes(width, "big") + data.b.to_bytes(width, "big")
    # điểm đường cong (curve.Point): nhận diện theo thuộc tính
    if hasattr(data, "is_infinity") and hasattr(data, "E"):
        if data.is_infinity():
            return b"P0"
        return b"P1" + encode(data.x) + encode(data.y)
    if isinstance(data, (list, tuple)):
        out = b"L" + _u32(len(data))
        for item in data:
            enc = encode(item)
            out += _u32(len(enc)) + enc
        return out
    raise TypeError("không mã hóa được kiểu %s" % type(data).__name__)


def _as_label(label):
    return label.encode("utf-8") if isinstance(label, str) else bytes(label)


class Transcript:
    """Transcript Fiat–Shamir; cùng chuỗi thao tác cho cùng thách thức."""

    def __init__(self, label, r=449):
        self.F = GF(r)
        lab = _as_label(label)
        self.state = hashlib.sha256(DOMAIN + _u32(len(lab)) + lab).digest()

    def _update(self, tag, label, payload):
        lab = _as_label(label)
        h = hashlib.sha256()
        h.update(self.state + tag + _u32(len(lab)) + lab + _u32(len(payload)) + payload)
        self.state = h.digest()

    def absorb(self, label, data):
        """Hấp thụ dữ liệu (bytes, int, phần tử trường, điểm, list...) kèm nhãn."""
        self._update(b"A", label, encode(data))

    def challenge_bytes(self, label, n=64):
        """Rút n byte giả ngẫu nhiên (n ≤ 64 ở mức này là đủ)."""
        self._update(b"C", label, b"")
        out, ctr = b"", 0
        while len(out) < n:
            out += hashlib.sha256(self.state + b"C" + bytes([ctr])).digest()
            ctr += 1
        # cập nhật trạng thái để thách thức sau khác thách thức trước
        self._update(b"D", label, out[:n])
        return out[:n]

    def challenge(self, label):
        """Thách thức trong F_r: 64 byte băm, đọc big-endian, lấy mod r."""
        return self.F(int.from_bytes(self.challenge_bytes(label, 64), "big"))
