function T = G2T_stiefel(G)
%G2R_STIEFEL Given G \in real^{nrs+1, d+1, N}, extract the N (nrs \times 1)
% vectors in the last column of the matrices spanning along the 3rd
% dimension

nrs = size(G,1)-1;
d = size(G,2)-1;
N = size(G,3);
T = zeros(nrs, d, N);

for ii = 1:size(G, 3)
    T(:,ii) = G(1:nrs, d+1, ii);
end


end

