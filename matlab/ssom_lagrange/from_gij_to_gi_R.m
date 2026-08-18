function [giR, booleans_R] = from_gij_to_gi_R(gijR, edges, N, R1_offset)

nrs = size(gijR, 1);
d = size(gijR, 2);
booleans_R = boolean(0) * ones(N,1); % alg should stop when all these are 1
booleans_R(1) = boolean(1); % node 1 chosen as reference
% d = size(T_diffs, 1);
giR = eye3d(nrs, d, N);
adjmat = edges2adjmatrix(edges);
g = digraph(adjmat);
giR(:,:,1) = R1_offset;
for ii = 2:N
    [shortest_p, length, edge_path] = shortestpath(g, ii, 1);
%     fprintf("shortest_p for node %g\n", ii);
%     disp(shortest_p)
%     fprintf("length for node %g\n", ii);
%     disp(length)
%     fprintf("edge_path for node %g\n", ii);
%     disp(edge_path)
    giR(:,:,ii) = R1_offset;
    for ep = flip(edge_path)
        giR(:,:,ii) = giR(:,:,ii) * gijR(:,:,ep);
    end
    % giR(:,:,ii) = giR(:,:,ii);
end
    

end