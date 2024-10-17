function adj_mat = make_adj_mat_from_edges(e, num_nodes)
%MAKE_ADJ_MAT_FROM_EDGES: Returns adjacency matrix computed from an edges
%matrix; the matrix has one edge per each row

adj_mat = zeros(num_nodes, num_nodes);
for kk = 1:size(e,1)
    ii = e(kk, 1);
    jj = e(kk, 2);
    adj_mat(ii, jj) = 1;
end

end %function
