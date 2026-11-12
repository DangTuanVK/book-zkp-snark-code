// Testbench chung cho montmul_* và ntt_bfly_*: đọc vector (sinh bởi
// ch11/sinh_vector_rtl.py từ snarklib), bơm một vector mỗi chu kỳ (chèn chu kỳ
// trống giả ngẫu nhiên để thử tín hiệu valid), so khớp đầu ra theo thứ tự.
// Cấu hình (W, số vector, file, tên DUT) lấy từ cfg.svh do rtl/run_all.sh chép vào.
`timescale 1ns/1ps
`include "cfg.svh"
module tb_rtl;
  localparam int W    = `TB_W;
  localparam int NVEC = `TB_NVEC;
`ifdef TB_MONTMUL
  localparam int NF = 3;                 // {a, b, y}
`else
  localparam int NF = 5;                 // {u, v, w_m, x, y}
`endif
  logic [NF*W-1:0] mem [0:NVEC-1];
  logic clk = 1'b0, rst_n = 1'b0, in_valid = 1'b0, out_valid;
  logic [W-1:0] i0, i1, i2, o0, o1;
  int sent = 0, got = 0, errors = 0, cycles = 0, first_in = -1, first_out = -1;
  int lfsr = 32'h1;

  always #5 clk = ~clk;

`ifdef TB_MONTMUL
  `TB_DUT dut (.clk, .rst_n, .in_valid, .a(i0), .b(i1), .out_valid, .y(o0));
  assign o1 = '0;
  assign i2 = '0;
`else
  `TB_DUT dut (.clk, .rst_n, .in_valid, .u(i0), .v(i1), .w_m(i2), .out_valid, .x(o0), .y(o1));
`endif

  // bơm vector: khoảng 1/8 số chu kỳ để trống (in_valid = 0)
  always @(posedge clk) begin
    if (rst_n) begin
      cycles <= cycles + 1;
      lfsr   <= {lfsr[30:0], lfsr[31] ^ lfsr[21] ^ lfsr[1] ^ lfsr[0]};
      if (sent < NVEC && lfsr[2:0] != 3'b000) begin
        in_valid <= 1'b1;
        i0 <= mem[sent][(NF-1)*W +: W];
        i1 <= mem[sent][(NF-2)*W +: W];
`ifndef TB_MONTMUL
        i2 <= mem[sent][(NF-3)*W +: W];
`endif
        sent <= sent + 1;
      end else begin
        in_valid <= 1'b0;
        i0 <= 'x; i1 <= 'x;             // dữ liệu rác khi không hợp lệ
`ifndef TB_MONTMUL
        i2 <= 'x;
`endif
      end
    end
  end

  // so khớp: đầu ra hợp lệ thứ k ứng với vector thứ k
  always @(posedge clk) begin
    if (rst_n && in_valid && first_in < 0) first_in = cycles;
    if (rst_n && out_valid) begin
      logic [W-1:0] e0, e1;
`ifdef TB_MONTMUL
      e0 = mem[got][0 +: W]; e1 = '0;
`else
      e0 = mem[got][W +: W]; e1 = mem[got][0 +: W];
`endif
      if (first_out < 0) first_out = cycles;
      if (o0 !== e0 || o1 !== e1) begin
        errors++;
        if (errors <= 5)
          $display("LOI vector %0d: ra (%h, %h), mong doi (%h, %h)", got, o0, o1, e0, e1);
      end
      got++;
    end
  end

  initial begin
    $readmemh(`TB_VEC, mem);
    repeat (3) @(posedge clk);
    rst_n = 1'b1;
    wait (got == NVEC);
    @(posedge clk);
    $display("%s: W = %0d, da kiem %0d vector, so loi = %0d, so chu ky = %0d, do tre = %0d chu ky",
             `TB_NAME, W, got, errors, cycles, first_out - first_in);
    if (errors == 0) $display("KET QUA: DAT");
    else             $display("KET QUA: TRUOT");
    $finish;
  end

  initial begin                           // chặn treo
    #(64'd10 * NVEC * 2 + 1000);
    $display("%s: HET GIO, moi nhan %0d/%0d vector", `TB_NAME, got, NVEC);
    $finish;
  end
endmodule
