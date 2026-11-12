// Testbench vét cạn: mọi cặp (a, b) trong [0,448]^2, so với (a+b) mod 449.
`timescale 1ns/1ps
module tb_addmod449;
  localparam int R = 449;
  logic [8:0] a, b, y;
  int errors = 0, total = 0;
  addmod449 dut (.a(a), .b(b), .y(y));
  initial begin
    for (int i = 0; i < R; i++) begin
      for (int j = 0; j < R; j++) begin
        a = i[8:0]; b = j[8:0];
        #1;
        total++;
        if (y !== (i + j) % R) begin
          errors++;
          if (errors <= 5) $display("LOI: a=%0d b=%0d y=%0d mong doi=%0d", i, j, y, (i+j)%R);
        end
      end
    end
    $display("addmod449: da kiem %0d cap, so loi = %0d", total, errors);
    $finish;
  end
endmodule
