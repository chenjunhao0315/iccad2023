module top(a0, a1, b1, b0, c, h0, h1, m0, m1);
input a0, a1, b1, b0, c;
output h0, h1, m0, m1;
wire a0t, b1t, n1, n2, n3, m0t;
and (h0, a0, b0);
not (b1t, b1);
not (a0t, a0);
nand (n1, a1, b0);
or (n2, b1t, a0t);
xor (h1, n1, n2);
nand (n3, b0, c);
xor (m0t, b1, n3);
not (m0, m0t);
xor (m1, b0, c);
endmodule
