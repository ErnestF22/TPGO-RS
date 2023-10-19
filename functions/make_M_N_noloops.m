function [Mmat_nl, Nmat_nl] = make_M_N_noloops(T_globalframe, Tijs_vec, edges, N)
%MAKE_M_N function that return M, N matrices, computing them without using
%any for loops; inputs can be noisy

num_edges = size(edges, 1);
d = size(T_globalframe, 1);
N = size(T_globalframe, 2);

%'nl' stands for 'no loops'

% Mmat_nl = zeros(d, num_edges);
% Nmat_nl = zeros(d, num_edges);

Mmat_nl = Tijs_vec;

Ti_pairs_mat = repmat(T_globalframe(:), 1, N) - repmat(T_globalframe, N, 1);

Ti_pairs_mat_cells = mat2cell(Ti_pairs_mat, d.*ones(N,1), ones(1,N));
Ti_pairs_mat_cells_vec = Ti_pairs_mat_cells(:);

edges_x = edges(:,1);
edges_y = edges(:,2);
edges_tot_id = N.*(edges_y - 1) + edges_x;

Nmat_nl_cells_vert = Ti_pairs_mat_cells_vec(edges_tot_id);
Nmat_nl_cells = reshape(Nmat_nl_cells_vert, 1, num_edges);

Nmat_nl = cell2mat(Nmat_nl_cells);

end %function

