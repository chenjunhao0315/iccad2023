module top( a, b, c, h1);
input a, b, c;
output h1;
wire n1, n2;
and (n1, b, c);
not (n2, n1);
and (h1, a, n2);
endmodule
