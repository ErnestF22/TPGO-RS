function [Mmat, Nmat] = make_M_N(T_globalframe, Tijs_vec, edges, ~)
%MAKE_M_N function that computes M, N matrices.
%Inputs can be noisy

num_edges = size(edges, 1);
d = size(T_globalframe, 1);

Mmat = zeros(d, num_edges);
Nmat = zeros(d, num_edges);

for edge_id = 1:num_edges
    ii = edges(edge_id, 1);
    jj = edges(edge_id, 2);
    Mmat(:,edge_id) = Tijs_vec(:, edge_id);
    Nmat(:,edge_id) = T_globalframe(:,ii) - T_globalframe(:,jj);
end

end %function