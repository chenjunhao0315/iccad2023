module top(x0, x1, y1, y0, z, u0, u1, w0, w1);
input x0, x1, y1, y0, z;
output u0, u1, w0, w1;
wire n1, n2, n3, w0t;
and (u0, x0, y0);
nand (n1, x1, y0);
nand (n2, y1, x0);
xor (u1, n1, n2);
nand (n3, y0, z);
xor (w0t, y1, n3);
not (w0, w0t);
xor (w1, y0, z);
endmodule
