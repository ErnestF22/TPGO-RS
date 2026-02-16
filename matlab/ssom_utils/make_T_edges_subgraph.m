function [T_edges] = make_T_edges_subgraph(T, edges, subgraph)
% MAKE_T_EDGES Return T_edges array with differences T(:,i) - T(:,j)
% for all (i,j) edge pairs of 'edges' input array


num_edges = size(edges, 1);
nrs = size(T, 1);
% N = size(T, 2);

% T_edges = zeros(nrs, num_edges);
id_edge_subgraph = 1;
for e = 1:num_edges
    
    ii = edges(e,1);
    jj = edges(e,2);
    if ismember(ii, subgraph) && ismember(jj, subgraph)
        T_edges(:,id_edge_subgraph) = T(:,ii) - T(:,jj); %TODO: fix preallocation
        id_edge_subgraph = id_edge_subgraph + 1; 
    end
end

% T1_offset = T(:,1);

end %file function
