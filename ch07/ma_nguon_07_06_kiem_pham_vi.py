"""Chương 7 -- Kiểm phạm vi theo cận tùy ý bằng LogUp (api cho Ch.12).

(1) logup.range_check(values, B): đồng nhất thức LogUp với bảng {0..B}.
(2) logup.vet_can(N, bounds): mọi (a, b) ∈ F_449^2 thỏa a·b = N, a ≤ Ba, b ≤ Bb.
(3) logup.prove_factor: proof mini-PLONK + LogUp cho từng nghiệm tìm được.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib import logup

# (1) cận tùy ý: -1 ≡ 448 và 20 nằm ngoài {0..19}
for vals, B in [([17, 19, 0], 19), ([17, 20], 19), ([17, 448], 19), ([23, 23, 5], 23),
                ([30, 28], 31), ([30, 28], 23)]:
    ok, info = logup.range_check(vals, B)
    print("range_check(%-12s, %2d) = %-5s ngoài bảng: %s" % (vals, B, ok, info["outside"]))
    assert ok == (not info["outside"])

# (2) vét cạn nghiệm của mạch có lookup
nghiem = {}
for bounds in [(31, 31), (19, 23)]:
    nghiem[bounds] = logup.vet_can(391, bounds)
    print("bounds = %s: %d nghiệm %s" % (bounds, len(nghiem[bounds]), nghiem[bounds]))

# (3) proof cho từng nghiệm; nghiệm của bảng (31, 31) thử lại với bảng (19, 23)
for a, b in nghiem[(31, 31)]:
    ok31 = logup.prove_factor(391, a, b, (31, 31), rng=random.Random(3))[0]
    try:
        ok_sua = logup.prove_factor(391, a, b, (19, 23), rng=random.Random(3))[0]
    except ValueError:
        ok_sua = "Prover dừng (ngoài bảng)"
    print("  (%2d, %2d): a·b = %3d; bảng (31, 31) -> %s; bảng (19, 23) -> %s"
          % (a, b, a * b, ok31, ok_sua))
    assert ok31 and (ok_sua is True) == ((a, b) == (17, 23))
    assert logup.accepted(391, a, b, (19, 23)) == ((a, b) == (17, 23))

assert nghiem[(31, 31)] == [(17, 23), (23, 17), (28, 30), (30, 28)]
assert nghiem[(19, 23)] == [(17, 23)] and 19 * 23 == 437 < 449
