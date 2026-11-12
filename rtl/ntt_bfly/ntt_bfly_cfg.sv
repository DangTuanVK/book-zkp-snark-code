// sinh bởi ch11/sinh_vector_rtl.py -- không sửa tay
`timescale 1ns/1ps
module ntt_bfly_f449 (input logic clk, rst_n, in_valid,
  input logic [15:0] u, v, w_m, output logic out_valid, output logic [15:0] x, y);
  ntt_bfly #(.W(16), .P(16'h1c1), .PINV(16'hf1bf)) u_bf (.*);
endmodule
module ntt_bfly_gl64 (input logic clk, rst_n, in_valid,
  input logic [63:0] u, v, w_m, output logic out_valid, output logic [63:0] x, y);
  ntt_bfly #(.W(64), .P(64'hffffffff00000001), .PINV(64'hfffffffeffffffff)) u_bf (.*);
endmodule
