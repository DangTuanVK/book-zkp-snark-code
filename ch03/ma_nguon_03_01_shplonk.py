"""KZG một điểm, batch một điểm và SHPLONK trên cùng bộ đa thức (lst:ch3-shplonk).

Bốn đa thức bậc 15 trên F_449; ba đa thức đầu mở tại z0 = 7, đa thức cuối
mở tại cả z0 và z1 = z0·ω_16 (giống z(X) mở tại zω trong PLONK). Prover và
verifier SHPLONK nằm trong snarklib/kzg.py (open_shplonk, verify_shplonk).
Đường cong toy là Type 1 (G2 = G1); τ = 77 chỉ để minh họa.
"""
import os
import random
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from snarklib import kzg, opcount as oc
from snarklib.poly import Poly

F = kzg.F
rng = random.Random(3)
srs = kzg.setup(15, 77)
fs = [Poly(F, [rng.randrange(449) for _ in range(16)]) for _ in range(4)]
Cs = [kzg.commit(srs, f) for f in fs]
z0 = F(7)
z1 = z0 * F.root_of_unity(16)                 # 7·349 mod 449
Ss = [[z0], [z0], [z0], [z0, z1]]

# (1) KZG một điểm cho đa thức thứ nhất
y, W = kzg.open(srs, fs[0], z0)
print("KZG một điểm: f_1(7) = %d, W = %s, verify = %s" % (int(y), W, kzg.verify(srs, Cs[0], z0, y, W)))

# (2) batch ba đa thức đầu tại z0: một nhân chứng
ys, Wb = kzg.open_batch(srs, fs[:3], Cs[:3], z0, kzg.fresh_transcript(b"vd-batch"))
okb = kzg.verify_batch(srs, Cs[:3], z0, ys, Wb, kzg.fresh_transcript(b"vd-batch"))
gb = kzg._gamma_one_point(kzg.fresh_transcript(b"vd-batch"), Cs[:3], z0, ys)
print("Batch 1 điểm: gamma = %d, y =" % int(gb), [int(v) for v in ys], " W =", Wb, " verify =", okb)

# (3) SHPLONK cho cả bốn đa thức, hai điểm: proof chỉ gồm (W, W')
with oc.scope() as cp:
    pr = kzg.open_shplonk(srs, fs, Cs, Ss, kzg.fresh_transcript(b"vd-shplonk"))
with oc.scope() as cv:
    ok = kzg.verify_shplonk(srs, Cs, Ss, pr, kzg.fresh_transcript(b"vd-shplonk"))
# tái tạo hai thách thức của verifier (γ trước W, z sau W) để in ra
tv = kzg.fresh_transcript(b"vd-shplonk")
kzg._absorb_claims(tv, Cs, Ss, pr.ys)
gs = tv.challenge(b"gamma")
tv.absorb(b"W", pr.W)
zs = kzg._shplonk_challenge_z(tv, [z0, z1])
print("SHPLONK: z1 = %d, giá trị mở =" % int(z1), [[int(v) for v in yy] for yy in pr.ys])
print("         gamma = %d, z = %d, Z_T(z) = %d" % (int(gs), int(zs), int(kzg.vanishing([z0, z1])(zs))))
print("         W = %s, W' = %s, verify = %s" % (pr.W, pr.Wp, ok))
print("         prover: MSM cỡ", cp["msm_sizes"], "; verifier:", oc.summary(cv))

# (4) gian lận: khai sai f_4(z1) -> bị từ chối
gian = pr._replace(ys=pr.ys[:3] + [[pr.ys[3][0], pr.ys[3][1] + 1]])
ok_gian = kzg.verify_shplonk(srs, Cs, Ss, gian, kzg.fresh_transcript(b"vd-shplonk"))
print("SHPLONK với f_4(z1) khai sai: verify =", ok_gian)
assert okb and ok and not ok_gian and cv["pairing"] == 2
