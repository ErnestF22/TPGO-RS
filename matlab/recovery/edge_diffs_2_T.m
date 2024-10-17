function [T, booleans_T] = edge_diffs_2_T(T_diffs, edges, N, T1_offset)
% MAKE_T_EDGES Return T_edges array with differences T(:,i) - T(:,j)
% for all (i,j) edge pairs of 'edges' input array
%

% num_edges = size(edges, 1);
nrs = size(T_diffs, 1);
booleans_T = boolean(0) * ones(N,1); % alg should stop when all these are 1
booleans_T(1) = boolean(1); % node 1 chosen as reference
% d = size(T_diffs, 1);
T = zeros(nrs, N);
adjmat = edges2adjmatrix(edges);
g = digraph(adjmat);
for ii = 2:N
    [shortest_p, length, edge_path] = shortestpath(g, ii, 1);
%     fprintf("shortest_p for node %g\n", ii);
%     disp(shortest_p)
%     fprintf("length for node %g\n", ii);
%     disp(length)
%     fprintf("edge_path for node %g\n", ii);
%     disp(edge_path)
    for ep = edge_path
        T(:,ii) = T(:,ii) + T_diffs(:,ep);
        disp('');
    end
    T(:,ii) = -T(:,ii);    
end

end %file function
