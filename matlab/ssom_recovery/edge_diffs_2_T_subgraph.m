function [T] = edge_diffs_2_T_subgraph(T_diffs_subgraph, edges, subgraph, offset)
% MAKE_T_EDGES Return T_edges array with differences T(:,i) - T(:,j)
% for all (i,j) edge pairs of 'edges' input array
%

% num_edges = size(edges, 1);
nrs = size(T_diffs_subgraph, 1);
% booleans_T = boolean(0) * ones(N,1); % alg should stop when all these are 1
% booleans_T(1) = boolean(1); % node 1 chosen as reference
% d = size(T_diffs, 1);
num_nodes_high_deg = size(subgraph, 2);
T = zeros(nrs, num_nodes_high_deg);

edges_in_subgraph_booleans0 = ismember(edges, subgraph);
edges_in_subgraph_booleans = ismember(edges_in_subgraph_booleans0, [1 1], "rows");

edges_subgraph = edges(edges_in_subgraph_booleans, :);
edges_subgraph_starting_from_1 = zeros(size(edges_subgraph));

for ii = 1:num_nodes_high_deg
    indices = edges_subgraph == subgraph(ii);
    edges_subgraph_starting_from_1(indices) = ii;
end

adjmat = edges2adjmatrix(edges_subgraph_starting_from_1); % !!
g = digraph(adjmat);

for ii = 2:num_nodes_high_deg
    [~, ~, edge_path] = shortestpath(g, ii, 1);
%     fprintf("shortest_p for node %g\n", ii);
%     disp(shortest_p)
%     fprintf("length for node %g\n", ii);
%     disp(length)
%     fprintf("edge_path for node %g\n", ii);
%     disp(edge_path)
    for ep = edge_path
        % ep_real = subgraph(1, ep);
        T(:,ii) = T(:,ii) + T_diffs_subgraph(:,ep);
        disp('');
    end
    T(:,ii) = -T(:,ii);    
end

T = T+offset;

end %file function
