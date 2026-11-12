#!/bin/bash
# Chạy lại toàn bộ phần RTL của Chương 11:
#   1) sinh vector + wrapper từ snarklib   (ch11/sinh_vector_rtl.py)
#   2) mô phỏng VCS mọi cấu hình            (rtl/tb/run_sim.sh -> ket_qua/ch11/sim_*.log)
#   3) tổng hợp DC lsi_10k có quét chu kỳ  (rtl/dc/run_dc.sh)
#   4) gom bảng                             (ch11/tong_hop_dc.py -> ket_qua/ch11/tong_hop.csv)
# Cách dùng: rtl/run_all.sh            (đủ, montmul_w256 tổng hợp lâu nhất)
#            rtl/run_all.sh --no-dc    (chỉ sinh vector và mô phỏng, khoảng 10 s)
set -e
RTL=$(cd "$(dirname "$0")" && pwd)
GOC=$(dirname "$RTL")
. "$RTL/env.sh"
mkdir -p "$GOC/ket_qua/ch11"
"$PY" "$GOC/ch11/sinh_vector_rtl.py" | tee "$GOC/ket_qua/ch11/sinh_vector.txt"
"$RTL/tb/run_sim.sh" | tee "$GOC/ket_qua/ch11/mo_phong.txt"
[ "$1" = "--no-dc" ] && exit 0
"$RTL/dc/run_dc.sh" | tee "$GOC/ket_qua/ch11/dc_quet.txt"
"$PY" "$GOC/ch11/tong_hop_dc.py" | tee "$GOC/ket_qua/ch11/tong_hop.txt"
