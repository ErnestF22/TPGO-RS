% this function compute A_cumulative and b_cumulative of all combinations
% of z:

function [A_c,b_c,z_set] = compute_A_cum_b_cum_all_combinations(Ab_set,L,n) 
% generate all possible z
set_all = table;
z_set=dec2bin(0:2^n-1)-'0';
set_all.z = z_set;

m = size(z_set,1);
A_c=[];
b_c=[];
for i=1:m
    z = z_set(i,:);
    [A,b] = get_cumulative(Ab_set,L,z);
    A_c = [A_c; A];
    b_c = [b_c; b];
end
Ab = [A_c,b_c];
idx = find(sum(abs(Ab),2)==0);
A_c(idx,:) = [] ;
b_c(idx,:) = [] ;
z_set(idx,:) = [] ;
end