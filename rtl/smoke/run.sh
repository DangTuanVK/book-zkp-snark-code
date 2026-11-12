#!/bin/bash
# Chạy mô phỏng VCS và tổng hợp DC cho addmod449. Mọi file sinh ra nằm trong out/.
set -e
cd "$(dirname "$0")"
export SNPSLMD_LICENSE_FILE=${SNPSLMD_LICENSE_FILE:-<may-chu-license>}
export VCS_HOME=${VCS_HOME:-<SYNOPSYS>/vcs/current}
VCS=<SYNOPSYS>/vcs/current/bin/vcs
DC=<SYNOPSYS>/syn/V-2023.12-SP3/bin/dc_shell
mkdir -p out
( cd out && $VCS -full64 -sverilog ../addmod449.sv ../tb_addmod449.sv -o simv -l vcs_compile.log >/dev/null && ./simv -l sim.log )
( cd out && /usr/bin/time -f "dc_shell wall=%e s" -o dc_time.txt $DC -f ../syn_addmod449.tcl > dc.log 2>&1 )
cat out/dc_time.txt
