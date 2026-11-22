"""Ch.11: gom kết quả DC (rtl/dc/out/<top>/) thành ket_qua/ch11/tong_hop.csv.

Mỗi lần compile (một chu kỳ ràng buộc) là một dòng. Cột chon = 1 đánh dấu lần
có chu kỳ nhỏ nhất mà vẫn đạt timing (slack >= 0) -- dòng dùng cho bảng của sách.
  tre_toi_han_ns = chu kỳ - slack (độ trễ đường tới hạn kể cả setup của flip-flop)
  f_uoc_luong_MHz = 1000 / tre_toi_han_ns
Báo cáo gốc (area, timing, qor) của dòng được chọn chép vào ket_qua/ch11/dc/.
Không phải mã in sách; rtl/run_all.sh gọi script này sau khi tổng hợp.
"""
import csv, os, re, shutil, sys

GOC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(GOC, "rtl", "dc", "out")
KQ = os.path.join(GOC, "ket_qua", "ch11")
THU_TU = ["montmul_w16", "montmul_w32", "montmul_w64", "montmul_w256",
          "ntt_bfly_f449", "ntt_bfly_gl64"]
MO_TA = {
    "montmul_w16": (16, "p = 449 (F_r của sách)"),
    "montmul_w32": (32, "BabyBear p = 15*2^27 + 1"),
    "montmul_w64": (64, "Goldilocks p = 2^64 - 2^32 + 1"),
    "montmul_w256": (256, "r của BLS12-381 (255 bit)"),
    "ntt_bfly_f449": (16, "butterfly trên F_449"),
    "ntt_bfly_gl64": (64, "butterfly trên Goldilocks"),
}


def so(rpt, nhan):
    m = re.search(re.escape(nhan) + r"\s*([0-9.]+)", rpt)
    return float(m.group(1)) if m else None


def main():
    os.makedirs(os.path.join(KQ, "dc"), exist_ok=True)
    dong = []
    for top in THU_TU:
        d = os.path.join(OUT, top)
        f = os.path.join(d, "ket_qua.txt")
        if not os.path.exists(f):
            print("bỏ qua %s (chưa tổng hợp)" % top)
            continue
        tt = open(os.path.join(d, "time.txt")).read().split()
        if not tt:
            print("bỏ qua %s (dc_shell chưa chạy xong)" % top)
            continue
        wall = tt[-1].replace("wall=", "")
        rows = []
        for ln in open(f):
            if ln.startswith("#") or not ln.strip():
                continue
            _, per, slack, arr, comp = ln.split()
            per, slack = float(per), float(slack)
            tag = "%s_T%s" % (top, ln.split()[1])
            area = open(os.path.join(d, tag + "_area.rpt")).read()
            qor = open(os.path.join(d, tag + "_qor.rpt")).read()
            tre = per - slack
            rows.append(dict(
                top=top, W=MO_TA[top][0], mo_ta=MO_TA[top][1], chu_ky_ns=per,
                slack_ns=round(slack, 3), dat=int(slack >= 0), tre_toi_han_ns=round(tre, 2),
                f_uoc_luong_MHz=round(1000 / tre, 2),
                dien_tich=so(area, "Total cell area:"), dt_to_hop=so(area, "Combinational area:"),
                dt_tuan_tu=so(area, "Noncombinational area:"),
                so_cell=int(so(area, "Number of cells:")),
                so_ff=int(so(area, "Number of sequential cells:")),
                muc_logic=int(so(qor, "Levels of Logic:")),
                compile_s=float(comp), dc_wall_tong_s=float(wall), chon=0, tag=tag))
        dat = [r for r in rows if r["dat"]]
        if dat:
            best = min(dat, key=lambda r: r["chu_ky_ns"])
            best["chon"] = 1
            for k in ("area", "timing", "qor"):
                shutil.copy(os.path.join(d, "%s_%s.rpt" % (best["tag"], k)),
                            os.path.join(KQ, "dc", "%s_%s.rpt" % (top, k)))
        dong += rows
    if not dong:
        # Bản clone mới không có rtl/dc/out/: giữ nguyên tong_hop.csv đã lưu.
        print("không có kết quả DC trong %s; giữ nguyên tong_hop.csv" % OUT)
        return
    cot = ["top", "W", "mo_ta", "chu_ky_ns", "slack_ns", "dat", "tre_toi_han_ns",
           "f_uoc_luong_MHz", "dien_tich", "dt_to_hop", "dt_tuan_tu", "so_cell", "so_ff",
           "muc_logic", "compile_s", "dc_wall_tong_s", "chon"]
    with open(os.path.join(KQ, "tong_hop.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cot, extrasaction="ignore")
        w.writeheader()
        w.writerows(dong)
    print("%-14s %4s %8s %8s %9s %10s %6s %5s %8s" % ("top", "W", "T(ns)", "tre(ns)", "f(MHz)",
                                                      "dien_tich", "cell", "FF", "dc(s)"))
    for r in dong:
        if r["chon"]:
            print("%-14s %4d %8.2f %8.2f %9.2f %10.0f %6d %5d %8.0f" % (
                r["top"], r["W"], r["chu_ky_ns"], r["tre_toi_han_ns"], r["f_uoc_luong_MHz"],
                r["dien_tich"], r["so_cell"], r["so_ff"], r["dc_wall_tong_s"]))


if __name__ == "__main__":
    main()
