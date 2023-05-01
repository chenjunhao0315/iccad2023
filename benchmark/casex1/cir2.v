module top( x, y, z, w, h2);
input x, y, z, w;
output h2;
wire n1, n2;
and (n1, x, y);
and (n2, z, w);
or (h2, n1, n2);
endmodule
