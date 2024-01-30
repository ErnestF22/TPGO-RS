function G_out = RT2G_stiefel(R,T)
%RT2G_STIEFEL Summary of this function goes here
%   Detailed explanation goes here

N=size(R,3);
d=size(R,2);
G_out=[R permute(T,[1 3 2])];
G_out=[G_out; zeros(1,d,N) ones(1,1,N)];


end

