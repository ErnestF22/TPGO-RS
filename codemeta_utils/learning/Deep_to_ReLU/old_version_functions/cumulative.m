function [layer,A_cum,b_cum]=cumulative(layer,L)
[m,n] = size(layer(1).A);
A_cum = eye(n);
b_cum = zeros(n,1);
for i=1:L
    Ai = layer(i).A;
    bi = layer(i).b;
    zi = layer(i).z;
    
    A_cum = zi.*Ai*A_cum;
    b_cum = zi.*Ai*b_cum+zi.*bi;
    
    layer(i).A_cum = A_cum;
    layer(i).b_cum = b_cum;
end
end