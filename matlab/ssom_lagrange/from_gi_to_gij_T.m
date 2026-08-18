function gijT = from_gi_to_gij_T(gi, edges)

nrs = 3; % non-Stiefel version
% N = size(giR, 3);

num_edges = size(edges, 1);

gijT = zeros(nrs, num_edges);

for ee = 1:num_edges
    ii = edges(ee, 1);
    jj = edges(ee, 2);
    tmp = inv(gi(:,:,ii)) * gi(:,:,jj);
    gijT(:,ee) = tmp(1:3, end);
end

end