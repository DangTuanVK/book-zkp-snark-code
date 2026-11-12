// Nhân Montgomery một từ, tham số hóa độ rộng W, R = 2^W.
//   Vào : a, b ở dạng Montgomery (a = ã, b = b̃), 0 <= a, b < P.
//   Ra  : y = REDC(a·b) = a·b·R^{-1} mod P.
// Kiến trúc: thanh ghi vào -> REDC tổ hợp trọn vẹn -> thanh ghi ra.
//   Độ trễ 2 chu kỳ, thông lượng 1 phép nhân/chu kỳ (nhận vector mới mỗi chu kỳ).
//   REDC tổ hợp gồm ba bộ nhân (a·b đủ 2W bit, m = T_lo·P' chỉ lấy W bit thấp,
//   m·P đủ 2W bit), một bộ cộng 2W+1 bit, phép dịch W bit (chỉ là chọn dây) và
//   một bước trừ có điều kiện.
// Lý do chọn: mỗi dòng RTL ứng đúng một dòng của REDC trong snarklib.montgomery,
//   dễ đối chiếu tay; phép chia cho R chỉ là chọn dây. Đổi lại đường tới hạn dài
//   (ba bộ nhân nối tiếp) -- chèn thêm tầng pipeline là bài tập cuối chương.
// Ràng buộc: P lẻ, P < 2^W, PINV = -P^{-1} mod 2^W. Khi a, b < P thì
//   T + m·P < 2·P·R nên t = (T + m·P) / R < 2P vừa W+1 bit; một lần trừ là đủ.
`timescale 1ns/1ps
module montmul #(
  parameter int W = 16,
  parameter logic [W-1:0] P    = 16'd449,
  parameter logic [W-1:0] PINV = 16'hf1bf     // -449^{-1} mod 2^16
) (
  input  logic         clk,
  input  logic         rst_n,
  input  logic         in_valid,
  input  logic [W-1:0] a,
  input  logic [W-1:0] b,
  output logic         out_valid,
  output logic [W-1:0] y
);
  logic [W-1:0]   a_q, b_q;     // thanh ghi vào
  logic           v_q;
  logic [2*W-1:0] T;            // T = ã·b̃
  logic [W-1:0]   m;            // m = (T mod R)·P' mod R
  logic [2*W:0]   U;            // U = T + m·P (chia hết cho R)
  logic [W:0]     t;            // t = U / R < 2P
  logic [W:0]     d;            // d = t - P
  logic [W-1:0]   y_c;

  always_ff @(posedge clk) begin
    a_q <= a;
    b_q <= b;
  end

  always_comb begin
    T   = a_q * b_q;                       // W x W -> 2W bit
    m   = T[W-1:0] * PINV;                 // ngữ cảnh W bit: tự cắt mod R
    U   = {1'b0, T} + m * P;               // ngữ cảnh 2W+1 bit
    t   = U[2*W:W];                        // chia cho R = lấy nửa trên
    d   = t - {1'b0, P};
    y_c = (t >= {1'b0, P}) ? d[W-1:0] : t[W-1:0];
  end

  always_ff @(posedge clk) y <= y_c;

  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      v_q       <= 1'b0;
      out_valid <= 1'b0;
    end else begin
      v_q       <= in_valid;
      out_valid <= v_q;
    end
  end
endmodule
