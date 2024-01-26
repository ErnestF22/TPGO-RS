function [Aconst, b] = make_Aconst_b(R_gf, Tijs_vec, edges)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nrs = size(R_gf, 1);
d = size(Tijs_vec, 1);
N = size(R_gf, 3);

adj_mat = make_adj_mat_from_edges(edges, N);

inc_mat = make_oriented_inc_mat_from_edges(edges, N);
inc_mat_rowsum = sum(inc_mat, 2);

% A

Aconst = repmat(inc_mat_rowsum', nrs, 1);

%b

% M = num2cell(adj_mat); %mask should be already implicit in Tijs_mat
Tijs_mat = tijs_vec_2_tijs_mat(Tijs_vec, edges, N);
Tijs_mat_cells = mat2cell(Tijs_mat, d*ones(1,N), ones(1,N));

% R_stacked = matStack(R_gf);
% R_cells = mat2cell(R_stacked, nrs*ones(1,N), d);
% rots_cells = repmat(R_cells, 1, N);

R_stacked = matStackH(R_gf);
R_cells = mat2cell(R_stacked, nrs, d*ones(1,N));
rots_cells = repmat(R_cells, N, 1);

b_cells = cellfun(@mtimes,rots_cells,Tijs_mat_cells,'Un',0); %apply mask
b_tmp = cell2mat(b_cells);
b = sum(b_tmp, 2);

end

