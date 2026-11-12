"""Chương 4 — PCS kiểu Ligero: commit và mở MLE của witness ví dụ N = 391 (lst:ch4-ligero).

Bảng 16 giá trị (quy ước z = (1, x, w) của snarklib, x = N, w = witness):
    z = (1, N, a, b, bit a_0..a_4, bit b_0..b_4, 0, 0),  a = 17, b = 23, bit LSB trước,
xem là MLE 4 biến. Ma trận 4 × 4, mã RS tỉ lệ 1/4 trên miền con 16 phần tử của F_449.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.field import GF
from snarklib.transcript import Transcript
from snarklib import ligero, mle

F = GF(449)
N, a, b = 391, 17, 23
bits = lambda v: [(v >> i) & 1 for i in range(5)]
z = [1, N, a, b] + bits(a) + bits(b) + [0, 0]
assert a * b == N and len(z) == 16
print("z =", z)

root, st = ligero.commit(z)                # mã hóa hàng + Merkle các cột
print("Ma trận M (4 x 4):", st["M"])
print("Hàng mã Enc(M[0]) (16 phần tử):", st["enc"][0])
print("Gốc Merkle:", root.hex()[:16], "...")

r = [F(2), F(3), F(5), F(7)]               # điểm mở (ở Chương 8 do sum-check sinh ra)
t = 4                                      # số cột mở
proof = ligero.open(st, r, Transcript(b"ch4-ligero"), t=t)
v = proof["v"]
print("r =", r, " v = z~(r) =", v, " (kiểm bằng mle.evaluate:", mle.evaluate(F, z, r), ")")
print("u_prox =", proof["u_prox"], " u_eval =", proof["u_eval"])
print("Cột được mở:", [j for j, _, _ in proof["cols"]])
for j, col, path in proof["cols"]:
    print("  cột", j, "=", col, " đường Merkle dài", len(path))

cnt = {}
ok = ligero.verify(root, r, v, proof, Transcript(b"ch4-ligero"), t=t, cnt=cnt)
print("Verifier chấp nhận:", ok, " chi phí:", cnt)
s = ligero.sizes(16, t=t)
print("Kích thước proof: %d phần tử trường + %d giá trị băm" % (s["proof_field"], s["proof_hash"]))

# Prover gian: sửa một phần tử trong cột được mở
j, col, path = proof["cols"][0]
forged = dict(proof, cols=[(j, [col[0] + 1] + col[1:], path)] + proof["cols"][1:])
bad = ligero.verify(root, r, v, forged, Transcript(b"ch4-ligero"), t=t)
print("Cột bị sửa -> chấp nhận:", bad)
assert ok and v == mle.evaluate(F, z, r) and not bad
