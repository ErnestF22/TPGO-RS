function gijR = from_gi_to_gij_R(giR, edges)
nrs = size(giR, 1);
d = size(giR, 2);
% N = size(giR, 3);

num_edges = size(edges, 1);

gijR = eye3d(nrs, d, num_edges);

for ee = 1:num_edges
    ii = edges(ee, 1);
    jj = edges(ee, 2);
    gijR(:,:,ee) = inv(giR(:,:,ii)) * giR(:,:,jj);
end

end
