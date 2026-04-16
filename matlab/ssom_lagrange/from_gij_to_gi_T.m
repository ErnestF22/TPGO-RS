function [giT, booleans_T] = from_gij_to_gi_T(gijT, gijR, edges, N, T1_offset, R1_offset)

nrs = size(gijR, 1);
d = size(gijR, 2);
booleans_T = boolean(0) * ones(N,1); % alg should stop when all these are 1
booleans_T(1) = boolean(1); % node 1 chosen as reference
% d = size(T_diffs, 1);
gi_full = zeros(nrs + 1, d + 1, N);
adjmat = edges2adjmatrix(edges);
g = digraph(adjmat);
gi_full(:,:,1) = RT2G(R1_offset, T1_offset);
for ii = 2:N
    [shortest_p, length, edge_path] = shortestpath(g, ii, 1);
%     fprintf("shortest_p for node %g\n", ii);
%     disp(shortest_p)
%     fprintf("length for node %g\n", ii);
%     disp(length)
%     fprintf("edge_path for node %g\n", ii);
%     disp(edge_path)
    gi_full(:,:,ii) = RT2G(R1_offset, T1_offset);
    for ep = flip(edge_path)
        tmp = RT2G(gijR(:,:,ep), gijT(:,ep));
        gi_full(:,:,ii) = gi_full(:,:,ii) * tmp;
    end
    % giR(:,:,ii) = giR(:,:,ii);
end

giT = G2T(gi_full);
    

end