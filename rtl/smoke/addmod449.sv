// Cộng modulo r = 449 (trường vô hướng F_r của sách), mạch tổ hợp.
// Giả thiết đầu vào đã rút gọn: a, b thuộc [0, 448].
`timescale 1ns/1ps
module addmod449 #(
  parameter int W = 9,
  parameter logic [W-1:0] R = 9'd449
) (
  input  logic [W-1:0] a,
  input  logic [W-1:0] b,
  output logic [W-1:0] y
);
  logic [W:0] s;   // tổng đủ 10 bit, tối đa 896
  logic [W:0] d;   // s - R (có thể âm)
  always_comb begin
    s = {1'b0, a} + {1'b0, b};
    d = s - {1'b0, R};
    // nếu s >= R thì lấy s - R, ngược lại giữ s
    y = (s >= {1'b0, R}) ? d[W-1:0] : s[W-1:0];
  end
endmodule
