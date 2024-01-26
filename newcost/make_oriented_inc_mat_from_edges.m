function incidence_mat = make_oriented_inc_mat_from_edges(e, num_nodes)
%MAKE_INDICENCE_MATRIX_FROM_EDGES: Returns indicence matrix computed from 
% an edges matrix; the matrix has one edge per each row
% 1 stands for start node, -1 for receiving node

    num_edges = size(e, 1);
    incidence_mat = zeros(num_nodes, num_edges);

    ctr = 1;
    for kk = 1:num_edges
        ii = e(kk, 1);
        jj = e(kk, 2);
        incidence_mat(ii, ctr) = 1;
        incidence_mat(jj, ctr) = -1;
        ctr = ctr + 1;
    end

end %function