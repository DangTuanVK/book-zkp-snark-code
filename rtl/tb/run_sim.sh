#!/bin/bash
# Mô phỏng VCS cho mọi cấu hình montmul và ntt_bfly; log thật chép vào ket_qua/ch11/.
# Cách dùng: rtl/tb/run_sim.sh [tên cấu hình ...]   (mặc định: tất cả)
# Cần chạy ch11/sinh_vector_rtl.py trước (rtl/run_all.sh làm việc này).
set -e
TB=$(cd "$(dirname "$0")" && pwd)
RTL=$(dirname "$TB")
KQ=$(dirname "$RTL")/ket_qua/ch11
. "$RTL/env.sh"
mkdir -p "$KQ"
CFGS=${*:-"montmul_w16 montmul_w32 montmul_w64 montmul_w256 ntt_bfly_f449 ntt_bfly_gl64"}
for c in $CFGS; do
  d=$TB/out/$c
  rm -rf "$d"; mkdir -p "$d"; cd "$d"
  cp "$TB/vec/cfg_$c.svh" cfg.svh
  ln -sf "$TB/vec/$c.hex" .
  $VCS -full64 -sverilog +incdir+. \
       "$RTL/montmul/montmul.sv" "$RTL/montmul/montmul_cfg.sv" \
       "$RTL/ntt_bfly/ntt_bfly.sv" "$RTL/ntt_bfly/ntt_bfly_cfg.sv" \
       "$TB/tb_rtl.sv" -o simv -l vcs_compile.log > /dev/null
  ./simv -l sim.log > /dev/null
  cp sim.log "$KQ/sim_$c.log"
  grep -E "da kiem|KET QUA|HET GIO" sim.log
  grep -q "KET QUA: DAT" sim.log || { echo "MO PHONG TRUOT: $c"; exit 1; }
done

# Kiểm tra testbench bắt được lỗi: bản montmul bỏ bước trừ cuối (đột biến), cấu hình W = 16.
# Kỳ vọng: TRUOT với đúng 155 lỗi = số cặp cần bước trừ cuối (ch11/ma_nguon_11_01_montmul_tay.py).
if [ -z "$*" ]; then
  d=$TB/out/dot_bien_montmul_w16
  rm -rf "$d"; mkdir -p "$d"; cd "$d"
  cp "$TB/vec/cfg_montmul_w16.svh" cfg.svh
  ln -sf "$TB/vec/montmul_w16.hex" .
  sed 's/y_c = (t >= {1.b0, P}) ? d\[W-1:0\] : t\[W-1:0\];/y_c = t[W-1:0];  \/\/ DOT BIEN/' \
      "$RTL/montmul/montmul.sv" > montmul_dot_bien.sv
  grep -q "DOT BIEN" montmul_dot_bien.sv
  $VCS -full64 -sverilog +incdir+. montmul_dot_bien.sv "$RTL/montmul/montmul_cfg.sv" \
       "$TB/tb_rtl.sv" -o simv -l vcs_compile.log > /dev/null
  ./simv -l sim.log > /dev/null || true
  cp sim.log "$KQ/sim_dot_bien_montmul_w16.log"
  echo "[dot bien: bo buoc tru cuoi]"; grep -E "da kiem|KET QUA" sim.log
  grep -q "so loi = 155," sim.log || { echo "DOT BIEN KHONG CHO 155 LOI"; exit 1; }
fi
