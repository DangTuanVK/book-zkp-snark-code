"""Đếm chi phí thật của các cách mở k đa thức bậc 15 trên F_449 (tab:ch3-chi-phi-batch).

Phương án (z0 = 7, z1 = z0·ω_16 như cặp z, zω của PLONK):
  rieng-1d   mở riêng từng đa thức tại z0          (k nhân chứng)
  batch-1d   batch k đa thức tại z0 bằng γ          (1 nhân chứng)
  rieng-2d   mở riêng từng đa thức tại z0 và z1     (2k nhân chứng, "ngây thơ")
  gwc-2d     mỗi điểm một nhân chứng, gộp pairing   (2 nhân chứng, PLONK §3.1)
  shplonk-2d SHPLONK cho cả hai điểm                (2 nhân chứng, BDFG20)
Số đếm lấy từ snarklib.opcount; không tính bước cam kết f_i (chung mọi phương án).
"""
import os
import random
import sys
import time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from snarklib import kzg, opcount as oc
from snarklib.poly import Poly

F = kzg.F
srs = kzg.setup(15, 77)
Z0 = F(7)
Z1 = Z0 * F.root_of_unity(16)


def chay(ten, fs, Cs):
    tr = lambda: kzg.fresh_transcript(b"dem-" + ten.encode())  # noqa: E731
    k = len(fs)
    with oc.scope() as cp:
        t0 = time.perf_counter()
        if ten == "rieng-1d":
            pf = [kzg.open(srs, f, Z0) for f in fs]
        elif ten == "batch-1d":
            pf = kzg.open_batch(srs, fs, Cs, Z0, tr())
        elif ten == "rieng-2d":
            pf = [kzg.open(srs, f, z) for f in fs for z in (Z0, Z1)]
        elif ten == "gwc-2d":
            pf = kzg.open_gwc(srs, fs, Cs, [[Z0, Z1]] * k, tr())
        else:
            pf = kzg.open_shplonk(srs, fs, Cs, [[Z0, Z1]] * k, tr())
        tp = time.perf_counter() - t0
    with oc.scope() as cv:
        t0 = time.perf_counter()
        if ten == "rieng-1d":
            ok = all(kzg.verify(srs, C, Z0, y, W) for C, (y, W) in zip(Cs, pf))
        elif ten == "batch-1d":
            ok = kzg.verify_batch(srs, Cs, Z0, pf[0], pf[1], tr())
        elif ten == "rieng-2d":
            ok = all(kzg.verify(srs, Cs[i // 2], (Z0, Z1)[i % 2], y, W) for i, (y, W) in enumerate(pf))
        elif ten == "gwc-2d":
            ok = kzg.verify_gwc(srs, Cs, [[Z0, Z1]] * k, pf, tr())
        else:
            ok = kzg.verify_shplonk(srs, Cs, [[Z0, Z1]] * k, pf, tr())
        tv = time.perf_counter() - t0
    assert ok, ten
    so_diem = 1 if ten.endswith("1d") else 2
    n_g1 = {"rieng-1d": k, "batch-1d": 1, "rieng-2d": 2 * k}.get(ten, 2)
    return dict(ten=ten, k=k, g1=n_g1, fr=k * so_diem, p_msm=cp["msm"],
                p_sizes=sorted(set(cp["msm_sizes"])), v_pair=cv["pairing"],
                v_fexp=cv["final_exp"], v_mul=cv.get("msm_nontrivial", 0) + cv.get("smul", 0),
                v_add=cv.get("add", 0), tp=tp * 1e3, tv=tv * 1e3)


if __name__ == "__main__":
    rng = random.Random(15)
    tieu_de = "%-11s %2s %4s %4s %6s %-9s %6s %5s %7s %8s %9s %9s" % (
        "phuong_an", "k", "G1", "F_r", "P_MSM", "co_MSM", "V_pair", "V_fe",
        "V_nhan", "V_cong", "t_P(ms)", "t_V(ms)")
    print(tieu_de)
    bang = []
    for k in (2, 4, 6, 8):
        fs = [Poly(F, [rng.randrange(449) for _ in range(16)]) for _ in range(k)]
        assert all(f.degree() == 15 for f in fs)
        Cs = [kzg.commit(srs, f) for f in fs]
        for ten in ("rieng-1d", "batch-1d", "rieng-2d", "gwc-2d", "shplonk-2d"):
            d = chay(ten, fs, Cs)
            bang.append(d)
            print("%-11s %2d %4d %4d %6d %-9s %6d %5d %7d %8d %9.1f %9.1f" % (
                d["ten"], k, d["g1"], d["fr"], d["p_msm"], ",".join(map(str, d["p_sizes"])),
                d["v_pair"], d["v_fexp"], d["v_mul"], d["v_add"], d["tp"], d["tv"]))
    # tự kiểm: số pairing đúng như công thức 2k / 2 / 4k / 2 / 2
    for d in bang:
        mong = {"rieng-1d": 2 * d["k"], "batch-1d": 2, "rieng-2d": 4 * d["k"]}.get(d["ten"], 2)
        assert d["v_pair"] == mong, d
    print("Chú thích: G1/F_r = số phần tử nhóm/trường trong proof; P_MSM = số MSM của prover "
          "(cỡ ghi ở co_MSM); V_pair = số pairing (vòng Miller), V_fe = số lũy thừa cuối; "
          "V_nhan = số phép nhân vô hướng G1 thật của verifier (vô hướng khác 0, 1); "
          "V_cong = số phép cộng/nhân đôi điểm của verifier; thời gian đo một lần, Python thuần.")
