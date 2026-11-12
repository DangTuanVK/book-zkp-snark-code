// Butterfly Cooley-Tukey (DIT) trên F_P: (x, y) = (u + ω·v, u - ω·v) mod P.
//   u, v ở dạng thường (hoặc cùng một dạng bất kỳ); twiddle w_m = ω·R mod P lưu
//   sẵn ở dạng Montgomery, nên REDC(v·w_m) = v·ω mod P ở đúng dạng của v --
//   không cần chuyển vào/ra dạng Montgomery trong vòng NTT.
// Pipeline 3 tầng, nhận một butterfly mới mỗi chu kỳ:
//   tầng 1-2 : montmul (thanh ghi vào, REDC tổ hợp, thanh ghi ra) -> t = ω·v
//              u đi song song qua hai thanh ghi trễ để gặp t đúng chu kỳ
//   tầng 3   : cộng và trừ modulo P song song, thanh ghi ra
`timescale 1ns/1ps
module ntt_bfly #(
  parameter int W = 16,
  parameter logic [W-1:0] P    = 16'd449,
  parameter logic [W-1:0] PINV = 16'hf1bf
) (
  input  logic         clk,
  input  logic         rst_n,
  input  logic         in_valid,
  input  logic [W-1:0] u,
  input  logic [W-1:0] v,
  input  logic [W-1:0] w_m,
  output logic         out_valid,
  output logic [W-1:0] x,
  output logic [W-1:0] y
);
  logic [W-1:0] t;              // t = ω·v mod P (ra khỏi montmul sau 2 chu kỳ)
  logic         t_valid;
  logic [W-1:0] u_d1, u_d2;     // u trễ 2 chu kỳ
  logic [W:0]   s, s_red;       // u + t, u + t - P
  logic [W:0]   dif, dif_cor;   // u - t, u - t + P
  logic [W-1:0] x_c, y_c;

  montmul #(.W(W), .P(P), .PINV(PINV)) u_mm (
    .clk(clk), .rst_n(rst_n), .in_valid(in_valid),
    .a(v), .b(w_m), .out_valid(t_valid), .y(t));

  always_ff @(posedge clk) begin
    u_d1 <= u;
    u_d2 <= u_d1;
  end

  always_comb begin
    s       = {1'b0, u_d2} + {1'b0, t};          // < 2P, vừa W+1 bit
    s_red   = s - {1'b0, P};
    x_c     = (s >= {1'b0, P}) ? s_red[W-1:0] : s[W-1:0];
    dif     = {1'b0, u_d2} - {1'b0, t};          // mượn nếu u < t
    dif_cor = dif + {1'b0, P};
    y_c     = (u_d2 >= t) ? dif[W-1:0] : dif_cor[W-1:0];
  end

  always_ff @(posedge clk) begin
    x <= x_c;
    y <= y_c;
  end

  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) out_valid <= 1'b0;
    else        out_valid <= t_valid;
  end
endmodule
