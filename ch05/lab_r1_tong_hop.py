"""Gom output thô của Lab R1 (ket_qua/ch05/lab_r1_raw.txt) thành lab_r1.csv và lab_r1.txt.

Chạy đo:  source rust/env.sh; RAYON_NUM_THREADS=16 $CARGO_TARGET_DIR/release/lab_r1_groth16_marlin 10,12,14,16 5
(và RAYON_NUM_THREADS=1); xem chú thích đầu rust/lab_r1_groth16_marlin/src/main.rs.
"""
import os

KQ = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "ket_qua", "ch05")
COT = ("he,k,ban_sao,rang_buoc,bien_instance,bien_witness,mien,so_lan,setup_ms,prove_ms,"
       "verify_ms,proof_byte,dung,tai_truoc,tai_sau").split(",")

rows, luong, ghi_chu = [], None, []
for line in open(os.path.join(KQ, "lab_r1_raw.txt"), encoding="utf-8"):
    if line.startswith("# Lab R1: RAYON_NUM_THREADS = "):
        luong = line.split("= ")[1].split(";")[0]
    elif line.startswith("# Lab R1 bắt đầu") or line.startswith("# rustc") or line.startswith("# kết thúc"):
        ghi_chu.append(line.strip())
    elif line.startswith("CSV,"):
        f = line.strip().split(",")[1:]
        d = dict(zip(COT, f[:len(COT)]))
        for extra in f[len(COT):]:
            k, v = extra.split("=")
            d[k] = v
        d["luong"] = luong
        rows.append(d)

head = ["he", "luong", "k", "ban_sao", "rang_buoc", "bien_instance", "bien_witness", "mien",
        "mien_K", "so_lan", "setup_ms", "prove_ms", "verify_ms", "srs_ms", "proof_byte",
        "dung", "tai_truoc", "tai_sau"]
with open(os.path.join(KQ, "lab_r1.csv"), "w", encoding="utf-8") as fo:
    fo.write(",".join(head) + "\n")
    for d in rows:
        fo.write(",".join(d.get(h, "") for h in head) + "\n")

out = ["# Lab R1: Groth16 (ark-groth16 0.6) và Marlin (ark-marlin 0.3), BLS12-381",
       "# Thời gian: trung vị của %s lần đo (ms). Máy: ket_qua/may_do.txt (máy dùng chung)." % rows[0]["so_lan"],
       "# setup: Groth16 = CRS riêng cho mạch; Marlin = index từ SRS phổ quát (SRS: cột srs, đo 1 lần)"]
out += ghi_chu
fmt = "%-8s %-5s %3s %6s %7s %7s %10s %10s %9s %9s %6s"
out.append(fmt % ("hệ", "luồng", "k", "bản sao", "ràng b.", "miền", "setup", "prove", "verify",
                  "srs", "byte"))
for d in rows:
    out.append(fmt % (d["he"], d["luong"], d["k"], d["ban_sao"], d["rang_buoc"], d["mien"],
                      d["setup_ms"], d["prove_ms"], d["verify_ms"], d.get("srs_ms", "-"),
                      d["proof_byte"]))
    assert d["dung"] == "true"
# tỉ lệ Marlin / Groth16 cùng số luồng, cùng k
idx = {(d["he"], d["luong"], d["k"]): d for d in rows}
out.append("# Tỉ lệ Marlin/Groth16 (prove), cùng k và số luồng:")
for (he, l, k), d in sorted(idx.items()):
    if he == "marlin" and ("groth16", l, k) in idx:
        g = idx[("groth16", l, k)]
        out.append("#   luồng %s, k = %s: prove x%.1f, setup/index x%.1f, verify x%.1f"
                   % (l, k, float(d["prove_ms"]) / float(g["prove_ms"]),
                      float(d["setup_ms"]) / float(g["setup_ms"]),
                      float(d["verify_ms"]) / float(g["verify_ms"])))
txt = "\n".join(out) + "\n"
open(os.path.join(KQ, "lab_r1.txt"), "w", encoding="utf-8").write(txt)
print(txt)
