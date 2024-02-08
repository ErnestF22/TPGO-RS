function G_out = RT2G_stiefel(R,T)
%RT2G_STIEFEL Compose affine "Stiefel-transformation" matrices

nrs=size(R,1);
d=size(R,2);
N=size(R,3);
G_out = zeros(nrs+1, d+1, N);


for ii = 1:N
    G_out(1:nrs, 1:d, ii) = R(:,:,ii);
    G_out(1:nrs, d+1, ii) = T(:,ii);
    G_out(nrs+1, d+1,ii) = 1; %"affine"
end



end

