KẾT QUẢ CHƯƠNG 11 (RTL: montmul, ntt_bfly) -- tái hiện: rtl/run_all.sh
Công cụ: VCS U-2023.03-SP2-2 (-full64 -sverilog), Design Compiler V-2023.12-SP3
(compile_ultra -no_autoungroup), thư viện GIÁO DỤC lsi_10k (không phải công nghệ thật:
đơn vị diện tích và ns của lsi_10k chỉ dùng để so sánh tương đối giữa các cấu hình).

sinh_vector.txt          tóm tắt vector sinh bằng snarklib (ch11/sinh_vector_rtl.py)
mo_phong.txt, sim_*.log  log VCS thật: số vector, số lỗi (= 0), độ trễ pipeline
ma_nguon_11_01_*.txt     mô hình bit-chính-xác montmul W = 16, ví dụ 17·23, số lần trừ cuối
ma_nguon_11_02_*.txt     mô hình butterfly F_449 và NTT n = 16 dựng từ butterfly
tong_hop.csv             mọi lần compile DC (quét chu kỳ); cột chon = 1 là chu kỳ nhỏ nhất đạt timing
tong_hop.txt             bảng rút gọn các dòng chon = 1
dc_quet.txt              output của rtl/dc/run_dc.sh (ket_qua.txt của từng wrapper + thời gian dc_shell)
dc/<top>_{area,timing,qor}.rpt   báo cáo DC gốc của dòng được chọn
dc/montmul_w256_T1000_long_*     lần thăm dò W = 256 với chu kỳ lỏng 1000 ns (chạy tay một lần,
                                 dc_shell 793 s), giữ lại để so diện tích theo ràng buộc

Định nghĩa cột tong_hop.csv:
  chu_ky_ns        chu kỳ ràng buộc (create_clock -period)
  slack_ns         slack của đường tới hạn (thanh ghi -> thanh ghi)
  tre_toi_han_ns   chu_ky_ns - slack_ns (gồm cả setup flip-flop 0,8 ns của lsi_10k)
  f_uoc_luong_MHz  1000 / tre_toi_han_ns -- ước lượng, không có wire load, không có layout
  dien_tich        Total cell area (đơn vị lsi_10k), dt_to_hop / dt_tuan_tu: tổ hợp / flip-flop
  muc_logic        Levels of Logic (report_qor)
  compile_s        thời gian riêng lệnh compile_ultra; dc_wall_tong_s: cả phiên dc_shell
                   (gồm khởi động, license, sinh alib ~55 s lần đầu và mọi lần quét)
