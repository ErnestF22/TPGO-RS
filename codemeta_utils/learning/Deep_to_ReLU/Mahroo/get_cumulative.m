function [A_cum,b_cum]=get_cumulative(Ab_set,L,z)
[~,n] = size(Ab_set(1).A);
A_cum = eye(n);
b_cum = zeros(n,1);
j=1;
for i=1:L
    Ai = Ab_set(i).A;
    bi = Ab_set(i).b;
    z_size = size(Ai,1);
    zi = z(:,j:j+z_size-1)';
    A_cum = zi.*Ai*A_cum;
    b_cum = zi.*Ai*b_cum+zi.*bi;
    j = j+z_size;
end
end