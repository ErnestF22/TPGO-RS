function m_i = get_Mi_from_Mmat(Mmat, edges, idx, d, N)
%GET_MI_FROM_MMAT Function that returns the M_i matrix used during
%Procrustes rotation estimation;
%in the variable, idx notation is used instead of i, because of the usual 
%problem with complex numbers and i, j in Matlab

m_i = zeros(d, N);
for edge_id = 1:size(edges, 1)
    if edges(edge_id, 1) == idx
        jj = edges(edge_id, 2);
        m_i(:, jj) = Mmat(:, edge_id);
    end

end

end %function