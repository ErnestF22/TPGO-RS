function R = G2R_stiefel(G)
%G2R_STIEFEL Given G \in real^{nrs+1, d+1, N}, extract the N (nrs \times d)
% submatrices spanning along the third dimension

nrs = size(G,1)-1;
d = size(G,2)-1;
N = size(G,3);
R = zeros(nrs, d, N);

for ii = 1:size(G, 3)
    R(:,:,ii) = G(1:nrs, 1:d, ii);
end


end

