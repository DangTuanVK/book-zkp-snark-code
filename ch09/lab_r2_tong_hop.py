"""Gom ket_qua/ch09/lab_r2.csv (Lab R2, Nova IVC) thành bảng cho tab:ch9-lab-r2 và fig:ch9-thoi-gian-buoc.

Chạy đo:  source rust/env.sh; cd rust/lab_r2_nova; RAYON_NUM_THREADS=16 cargo run --release
(output thật: ket_qua/ch09/lab_r2.txt, lab_r2.csv).
"""
import csv, os

KQ = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "ket_qua", "ch09")
rows = list(csv.DictReader(open(os.path.join(KQ, "lab_r2.csv"), encoding="utf-8")))
v = lambda x: ("%.1f" % float(x)).replace(".", ",")

print("Bảng 1 -- theo số bước (F = kiểm phân tích, 14 ràng buộc; mạch chính 9833, mạch phụ 10349)")
print("số bước | setup PP (ms) | bước TB / TV (ms) | tổng IVC (s) | verify RS (ms) | nén (ms) | verify nén (ms) | proof nén (byte)")
for r in rows:
    if r["them"] == "0":
        print("%7s | %s | %s / %s | %s | %s | %s | %s | %s" % (
            r["so_buoc"], v(r["setup_ms"]), v(r["buoc_tb_ms"]), v(r["buoc_tv_ms"]),
            v(float(r["tong_ivc_ms"]) / 1000), v(r["verify_rs_ms"]), v(r["nen_prove_ms"]),
            v(r["nen_verify_ms"]), r["nen_byte"]))
print("\nBảng 2 -- theo cỡ mạch bước (10 bước): ràng buộc thêm -> thời gian mỗi bước")
print("thêm | rb mạch chính | bước TV (ms) | nén (ms) | proof nén (byte) | RecursiveSNARK (byte)")
for r in rows:
    if r["so_buoc"] == "10":
        print("%6s | %6s | %s | %s | %s | %s" % (r["them"], r["rb_chinh"], v(r["buoc_tv_ms"]),
              v(r["nen_prove_ms"]), r["nen_byte"], r["rs_byte"]))
base = [r for r in rows if r["them"] == "0"]
assert len({r["nen_byte"] for r in base}) == 1, "proof nén phải không đổi theo số bước"
print("\nKích thước proof nén không đổi theo số bước:", base[0]["nen_byte"], "byte")
print("Tải máy (loadavg) trước/sau:", base[0]["tai_truoc"], "->", rows[-1]["tai_sau"])
