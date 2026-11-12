// sinh bởi ch11/sinh_vector_rtl.py -- không sửa tay
// Wrapper cố định p và p' = -p^{-1} mod 2^W cho từng cấu hình (dùng cho mô phỏng và DC).
`timescale 1ns/1ps
// montmul_w16: p = F_449 của sách, R = 2^16
module montmul_w16 (input logic clk, rst_n, in_valid,
  input logic [15:0] a, b, output logic out_valid, output logic [15:0] y);
  montmul #(.W(16), .P(16'h1c1), .PINV(16'hf1bf)) u_mm (.*);
endmodule
// montmul_w32: p = BabyBear 15·2^27 + 1 (31 bit)
module montmul_w32 (input logic clk, rst_n, in_valid,
  input logic [31:0] a, b, output logic out_valid, output logic [31:0] y);
  montmul #(.W(32), .P(32'h78000001), .PINV(32'h77ffffff)) u_mm (.*);
endmodule
// montmul_w64: p = Goldilocks 2^64 - 2^32 + 1
module montmul_w64 (input logic clk, rst_n, in_valid,
  input logic [63:0] a, b, output logic out_valid, output logic [63:0] y);
  montmul #(.W(64), .P(64'hffffffff00000001), .PINV(64'hfffffffeffffffff)) u_mm (.*);
endmodule
// montmul_w256: p = r của BLS12-381 (255 bit)
module montmul_w256 (input logic clk, rst_n, in_valid,
  input logic [255:0] a, b, output logic out_valid, output logic [255:0] y);
  montmul #(.W(256), .P(256'h73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001), .PINV(256'h3d443ab0d7bf2839181b2c170004ec0653ba5bfffffe5bfdfffffffeffffffff)) u_mm (.*);
endmodule
