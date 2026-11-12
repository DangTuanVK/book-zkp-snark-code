# Tổng hợp addmod449 bằng Design Compiler với thư viện giáo dục lsi_10k.
# Chạy từ thư mục out/ (run.sh làm việc này) để file phụ của DC không rơi vào mã nguồn.
if {[info exists ::env(SYN_LIB_DIR)] && $::env(SYN_LIB_DIR) ne ""} {
    set LIB_DIR $::env(SYN_LIB_DIR)
} else {
    error "Cần biến môi trường SYN_LIB_DIR (thư mục chứa lsi_10k.db); xem rtl/env.sh"
}
set search_path   [list . $LIB_DIR]
set target_library lsi_10k.db
set link_library   [list * lsi_10k.db]
analyze -format sverilog ../addmod449.sv
elaborate addmod449
current_design addmod449
link
# mạch tổ hợp: ràng buộc độ trễ tối đa từ vào tới ra (đơn vị ns của lsi_10k)
set_max_delay 10 -from [all_inputs] -to [all_outputs]
compile_ultra -no_autoungroup
report_area   > area.rpt
report_timing > timing.rpt
report_qor    > qor.rpt
write -format verilog -hierarchy -output addmod449_net.v
exit
