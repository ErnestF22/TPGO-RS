function [T_edges, T1_offset] = make_T_edges(T, edges)
% MAKE_T_EDGES Return T_edges array with differences T(:,i) - T(:,j)
% for all (i,j) edge pairs of 'edges' input array
%


num_edges = size(edges, 1);
nrs = size(T, 1);
% N = size(T, 2);
T_edges = zeros(nrs, num_edges);
for e = 1:num_edges
    ii = edges(e,1);
    jj = edges(e,2);
    T_edges(:,e) = T(:,ii) - T(:,jj);
end

T1_offset = T(:,1);

end %file function
