"""Chương 12 -- Fuzz witness: đột biến witness hợp lệ, đếm đột biến vẫn thỏa.

Hai mức: (1) R1CS -- sửa trực tiếp các ô của w; (2) mạch -- đột biến (a, b) và
sinh lại phần còn lại của witness bằng gợi ý. Seed cố định để output tái lập.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.dsl import mach_phan_tich, mach_phan_tich_sua
from snarklib import uniqueness as U

DUNG = {"a": 17, "b": 23}
for ten, c in (("mạch cũ", mach_phan_tich()[0]), ("mạch đã sửa", mach_phan_tich_sua()[0])):
    r1cs = c.compile_r1cs()
    x, w = c.r1cs_assignment(dict(N=391, **DUNG))
    print("== %s (m = %d) ==" % (ten, r1cs.m))
    # (1) mức R1CS: 5 kiểu đột biến x 2000 lần
    kq = U.fuzz_witness(r1cs, x, w, trials=2000, rng=random.Random(12))
    for kind, (n, hit, new) in kq.items():
        print("  R1CS %-15s: %5d lần, %d đột biến thỏa" % (kind, n, hit))
    # (2) mức mạch: (a, b) ngẫu nhiên trong [0, 31]^2, rồi trong cả F_449^2
    for mo_ta, pool in (("[0,31]^2", range(32)), ("F_449^2", range(449))):
        n, seen = U.fuzz_inputs(c, {"N": 391}, DUNG, trials=5000,
                                rng=random.Random(2024), spread={"a": pool, "b": pool})
        tong = sum(seen.values())
        print("  mạch (a,b) ∈ %-8s: %5d lần, %3d lần thỏa (%.2f%%), cặp mới: %s"
              % (mo_ta, n, tong, 100.0 * tong / n, sorted(seen)))
        if ten == "mạch cũ" and pool == range(32):
            assert sorted(seen) == [(23, 17), (28, 30), (30, 28)]
        if ten == "mạch đã sửa":
            assert not seen
    assert all(hit == 0 for _, hit, _ in kq.values())
print("Kết luận: đột biến từng ô R1CS không bao giờ lọt (các bit và tổ hợp ràng nhau);")
print("fuzz theo đầu vào tìm ra nghiệm giả khi miền thử nhỏ, gần như bỏ sót trên F_449^2.")
