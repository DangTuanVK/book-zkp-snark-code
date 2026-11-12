#!/bin/bash
# Tổng hợp DC (lsi_10k) cho các wrapper montmul_* và ntt_bfly_*, có quét chu kỳ.
# Cách dùng: rtl/dc/run_dc.sh [tên ...]   (mặc định: tất cả, kể cả montmul_w256)
# Mỗi wrapper chạy trong rtl/dc/out/<tên>/; thời gian dc_shell ghi vào time.txt.
# Gom kết quả: ch11/tong_hop_dc.py -> ket_qua/ch11/tong_hop.csv và báo cáo gốc.
set -e
DCDIR=$(cd "$(dirname "$0")" && pwd)
export SNPSLMD_LICENSE_FILE=${SNPSLMD_LICENSE_FILE:-<may-chu-license>}
DC=<SYNOPSYS>/syn/V-2023.12-SP3/bin/dc_shell

# tên -> "P0 NSWEEP NREFINE" (P0: chu kỳ khởi đầu, ns; chọn hơi lỏng so với độ trễ dự kiến)
declare -A CFG=(
  [montmul_w16]="50 6 2"
  [montmul_w32]="45 8 2"
  [montmul_w64]="90 8 2"
  [ntt_bfly_f449]="45 8 2"
  [ntt_bfly_gl64]="100 8 2"
  [montmul_w256]="120 1 0"
)
# montmul_w256: một lần compile ở 120 ns đã mất 3453 s (57,5 phút; dc_shell 3477 s), nên
# không quét; lần thăm dò chạy tay ở 1000 ns mất 793 s (ket_qua/ch11/dc/montmul_w256_T1000_long_*).
TOPS=${*:-"montmul_w16 montmul_w32 montmul_w64 ntt_bfly_f449 ntt_bfly_gl64 montmul_w256"}
for t in $TOPS; do
  read -r P0 NS NR <<< "${CFG[$t]}"
  d=$DCDIR/out/$t
  rm -rf "$d"; mkdir -p "$d"; cd "$d"
  /usr/bin/time -f "%e" -o time.txt \
    $DC -x "set TOP $t; set P0 $P0; set NSWEEP $NS; set NREFINE $NR" -f "$DCDIR/syn.tcl" > dc.log 2>&1
  echo "$t: dc_shell $(cat time.txt) s"; cat ket_qua.txt
  if grep -q "^Error" dc.log; then echo "DC BAO LOI: $t"; grep "^Error" dc.log; exit 1; fi
done
