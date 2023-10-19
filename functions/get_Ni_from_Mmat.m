function n_i = get_Ni_from_Mmat(Nmat, edges, idx, d, N)
%GET_MI_FROM_MMAT Function that returns the M_i matrix used during
%Procrustes rotation estimation;
%idx notation is used instead of i, because of the usual problem with
%complex numbers

n_i = zeros(d, N);
for edge_id = 1:size(edges, 1)
    if edges(edge_id, 1) == idx
        jj = edges(edge_id, 2);
        n_i(:, jj) = Nmat(:, edge_id);
    end

end

end %function