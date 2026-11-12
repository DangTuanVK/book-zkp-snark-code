"""Chương 4 — Đếm phép toán của Verifier IPA theo n = 2^3..2^8 (Mục 4.2a).

So sánh: Bulletproofs IPA tổng quát kiểm "gập từng vòng" (Protocol 2) với kiểm
"một MSM" (Bünz et al. Mục 3.1), và PCS từ IPA (b công khai, Halo).
smul = phép nhân vô hướng trên nhóm, add = phép cộng điểm, fmul = phép nhân trường
(chỉ tính các phép trong véc-tơ s và <s, b>).
"""
import csv, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from snarklib.transcript import Transcript
from snarklib.pedersen import PedersenParams, msm, padd
from snarklib import ipa

rows = []
for l in range(3, 9):
    n = 1 << l
    G, H = PedersenParams(n, b"dG").G, PedersenParams(n, b"dH").G
    U = PedersenParams(1, b"dU").G[0]
    a = [(5 * i + 1) % 449 for i in range(n)]
    b = [(3 * i * i + 2) % 449 for i in range(n)]
    P, c = padd(msm(a, G), msm(b, H)), ipa.inner(a, b)
    pr = ipa.bp_prove(G, H, U, a, b, Transcript(b"dem"))
    row = {"n": n}
    for m in ("fold", "msm"):
        k = ipa.VCount()
        assert ipa.bp_verify(G, H, U, P, c, pr, Transcript(b"dem"), m, k)
        row.update({"bp_%s_%s" % (m, key): val for key, val in k.as_dict().items()})
    pp = ipa.setup(n)
    C = ipa.commit(pp, a)
    v, pp_pr = ipa.open(pp, a, 7, Transcript(b"dem"))
    k = ipa.VCount()
    assert ipa.verify(pp, C, 7, v, pp_pr, Transcript(b"dem"), cnt=k)
    row.update({"pcs_" + key: val for key, val in k.as_dict().items()})
    rows.append(row)

print("%5s | %-15s | %-15s | %-15s" % ("n", "BP gập (smul)", "BP 1 MSM (smul)", "PCS (smul/fmul)"))
for r in rows:
    print("%5d | %-15d | %-15d | %d/%d" % (r["n"], r["bp_fold_smul"], r["bp_msm_smul"],
                                           r["pcs_smul"], r["pcs_fmul"]))
out = os.path.join(ROOT, "ket_qua", "ch04", "dem_phep_toan_ipa.csv")
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
print("Đã ghi", os.path.relpath(out, ROOT))
# Công thức: BP 1 MSM = 2 (x·U, (x·c)·U) + 2 log n (L, R) + 2n (<s,G>, <1/s,H>) + 3 (a·G, b·H, ab·U);
# BP gập = 4(n-1) + 2 log n + 5;  PCS = 2 (x_u·U, v·U') + 2 log n + n + 2
assert all(r["bp_msm_smul"] == 2 * r["n"] + 2 * (r["n"].bit_length() - 1) + 5 for r in rows)
assert all(r["bp_fold_smul"] == 4 * (r["n"] - 1) + 2 * (r["n"].bit_length() - 1) + 5 for r in rows)
assert all(r["pcs_smul"] == r["n"] + 2 * (r["n"].bit_length() - 1) + 4 for r in rows)
