function y = halfspaceout(x,A_out,b_out,A_cum,b_cum)
y=NaN;
if (A_out*x+b_out<=0)
    y=A_cum*x+b_cum;
end
end