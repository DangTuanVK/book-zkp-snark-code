#!/usr/bin/env bash
# Lab R3 (Ch.10): chạy lab_r3_msm_ntt hai lần (đa luồng mặc định, rồi 1 luồng),
# ghi ket_qua/ch10/r3.csv (phep,n,luong,trung_vi_ms,so_lan) và r3.txt (log đầy đủ).
# Cách chạy:  bash ch10/chay_r3.sh [log_min] [log_max]
set -euo pipefail
GOC="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$GOC/rust/env.sh"
LMIN="${1:-10}"; LMAX="${2:-20}"
OUT="$GOC/ket_qua/ch10"; mkdir -p "$OUT"
LOG="$OUT/r3.txt"; CSV="$OUT/r3.csv"
BIN="$CARGO_TARGET_DIR/release/lab_r3_msm_ntt"
(cd "$GOC/rust/lab_r3_msm_ntt" && cargo build --release --quiet)

{
  echo "# Lab R3 -- output thật, $(date '+%Y-%m-%d %H:%M:%S %z')"
  echo "# Máy: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //'), nproc = $(nproc)"
  echo "# $(rustc --version); arkworks 0.6.0 (ark-ec/ff/poly bật parallel)"
  echo "# Máy dùng chung: tải thay đổi giữa các lần đo, xem uptime trước/sau mỗi lượt."
} > "$LOG"

for MODE in multi single; do
  echo "" >> "$LOG"
  if [ "$MODE" = single ]; then
    echo "===== RAYON_NUM_THREADS=1 $BIN $LMIN $LMAX =====" >> "$LOG"
    export RAYON_NUM_THREADS=1
  else
    echo "===== (đa luồng mặc định) $BIN $LMIN $LMAX =====" >> "$LOG"
    unset RAYON_NUM_THREADS
  fi
  echo "uptime trước: $(uptime)" >> "$LOG"
  T0=$(date +%s)
  "$BIN" "$LMIN" "$LMAX" >> "$LOG" 2>&1
  echo "uptime sau:   $(uptime)  (lượt chạy $(( $(date +%s) - T0 )) s)" >> "$LOG"
done

echo "phep,n,luong,trung_vi_ms,so_lan" > "$CSV"
grep '^CSV,' "$LOG" | cut -d, -f2- >> "$CSV"
echo "Đã ghi $CSV ($(($(wc -l < "$CSV") - 1)) dòng) và $LOG"
