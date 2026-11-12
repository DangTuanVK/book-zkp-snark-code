"""Chương 4 — Kích thước proof IPA, Hyrax và Ligero theo n = 2^3..2^8 (fig:ch4-kich-thuoc-proof).

Mỗi n: chạy thật commit/open/verify trên bảng ngẫu nhiên tất định, ĐẾM phần tử
trong proof trả về. Quy đổi byte giả định 32 byte cho mỗi điểm nén, mỗi phần tử
trường và mỗi giá trị băm (cỡ đường cong/trường 255 bit, SHA-256) — trên toy
r = 449 các phần tử thật chỉ 2 byte. Ligero: ρ = 1/4, t = 8 cột mở (đối số dòng lệnh).
"""
import csv, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from snarklib.field import GF
from snarklib.transcript import Transcript
from snarklib import ipa, hyrax, ligero

F, B = GF(449), 32
T_COLS = int(sys.argv[1]) if len(sys.argv) > 1 else 8
rows = []
for l in range(3, 9):
    n = 1 << l
    tab = [(7 * i * i + 3 * i + 1) % 449 for i in range(n)]
    pt = [F(3 + 2 * i) for i in range(l)]
    # IPA (PCS đơn biến, b = (1, z, ..., z^{n-1}))
    pp = ipa.setup(n)
    C = ipa.commit(pp, tab)
    v, pr = ipa.open(pp, tab, 5, Transcript(b"kt"))
    c = ipa.VCount()
    assert ipa.verify(pp, C, 5, v, pr, Transcript(b"kt"), cnt=c)
    ig, iF = len(pr["L"]) + len(pr["R"]), 1
    # Hyrax (đa tuyến tính, ma trận √n × √n, mở bằng IPA)
    hp = hyrax.setup(n)
    cm, st = hyrax.commit(hp, tab)
    hv, hpr = hyrax.open(hp, st, pt, Transcript(b"kt"))
    assert hyrax.verify(hp, cm, pt, hv, hpr, Transcript(b"kt"), n)
    hg = len(hpr["ipa"]["L"]) * 2
    # Ligero (đa tuyến tính)
    root, ls = ligero.commit(tab)
    lp = ligero.open(ls, pt, Transcript(b"kt"), t=T_COLS)
    lc = {}
    assert ligero.verify(root, pt, lp["v"], lp, Transcript(b"kt"), t=T_COLS, cnt=lc)
    lf = len(lp["u_prox"]) + len(lp["u_eval"]) + sum(len(col) for _, col, _ in lp["cols"])
    lh = sum(len(path) for _, _, path in lp["cols"])
    rows.append(dict(n=n, log2n=l, ipa_group=ig, ipa_field=iF, ipa_bytes=B * (ig + iF),
                     ipa_ver_smul=c.smul, hyrax_commit_group=len(cm), hyrax_proof_group=hg,
                     hyrax_proof_field=1, hyrax_proof_bytes=B * (hg + 1),
                     ligero_m=ls["m"], ligero_k=ls["k"], ligero_t=len(lp["cols"]),
                     ligero_field=lf, ligero_hash=lh, ligero_bytes=B * (lf + lh),
                     ligero_ver_hash=lc["hash"], ligero_ver_fmul=lc["fmul"]))

cols = list(rows[0])
print("%5s %9s %8s %9s %9s %9s %10s %9s" % ("n", "IPA(G+F)", "IPA B", "IPA smul",
      "Hyrax B", "Lig(F+H)", "Ligero B", "Lig hash"))
for r in rows:
    print("%5d %5d+%-3d %8d %9d %9d %5d+%-4d %9d %9d" % (r["n"], r["ipa_group"],
          r["ipa_field"], r["ipa_bytes"], r["ipa_ver_smul"], r["hyrax_proof_bytes"],
          r["ligero_field"], r["ligero_hash"], r["ligero_bytes"], r["ligero_ver_hash"]))
out = os.path.join(ROOT, "ket_qua", "ch04", "kich_thuoc_proof.csv")
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    w.writerows(rows)
print("Đã ghi", os.path.relpath(out, ROOT))
assert all(r["ipa_group"] == 2 * r["log2n"] for r in rows)
assert all(r["ligero_bytes"] > r["ipa_bytes"] for r in rows)
