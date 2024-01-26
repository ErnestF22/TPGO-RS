function [LR, PR, BR_const] = make_LR_PR_BR_noloops(R_gf, Tijs, edges)
%MAKE_LR_PR_BR Summary of this function goes here
%   Detailed explanation goes here

nrs = size(R_gf, 1);
d = size(R_gf, 2);
N = size(R_gf, 3);

num_edges = size(edges, 1);

adj_mat = make_adj_mat_from_edges(edges, N);
M = num2cell(adj_mat); %mask

Tijs_mat = tijs_vec_2_tijs_mat(Tijs, edges, N);

%%
%LR
%LR_{ij} = b_{ij} * b_{ij}'
%LR = sum_{e \in (i,j)} LR_{ij}
B1A = zeros(N*N,N);
one_rows_indices = 1:N+1:N*N;
B1A(one_rows_indices, :) = ones(N,N);
B1B = - repmat(eye(N), N, 1);
B1 = B1A + B1B;

B1_cells = mat2cell(B1, N*ones(N,1), ones(1,N));
B1_cells_transp = cellfun(@transp, B1_cells, 'Un', 0);
BIJS = cellfun(@mtimes, B1_cells, B1_cells_transp, 'Un', 0);
BIJS_masked = cellfun(@mtimes, BIJS, M, 'Un', 0);
BIJS_masked_vectorized = cellfun(@vec, BIJS_masked, 'Un', 0);
%OBS. col-major order is not very important to be kept here
%since we are summing all cells
LR_vec = sum(reshape(cell2mat(BIJS_masked_vectorized'),N*N,[]), 2); %sum all
LR = reshape(LR_vec, N, N);

%%
%PR
%PR_{ij} = b_{ij} * T_{ij}' * R_i'
%PR = sum_{e \in (i,j)} PR_{ij}

Tijs_mat_cells = mat2cell(Tijs_mat, d*ones(N,1), ones(1,N));
Tijs_transp_mat_cells = cellfun(@transp, Tijs_mat_cells, 'Un', 0);

R_PR = repmat(matStack(multitransp(R_gf)), 1, N);
R_PR_cells = mat2cell(R_PR, d*ones(N,1), nrs*ones(1,N));

PR_12_cells = cellfun(@mtimes, B1_cells, Tijs_transp_mat_cells, 'Un', 0);
PR_123_cells = cellfun(@mtimes, PR_12_cells, R_PR_cells, 'Un', 0);
PR_123_cells_vec = cellfun(@vec, PR_123_cells, 'Un', 0);
PR_vec = sum(reshape(cell2mat(PR_123_cells_vec(:)),nrs*N,[]), 2); %sum all
PR = reshape(PR_vec, N, nrs);


%%
%BR
%BR_{ij} = T_{ij} * T_{ij}'
%BR = sum_{e \in (i,j)} BR_{ij}

BR_const_ijs = cellfun(@mtimes, ...
    Tijs_mat_cells, Tijs_transp_mat_cells, 'Un', 0);
BR_const_ijs_vec = cellfun(@vec, BR_const_ijs, 'Un', 0);
BR_const_vec = sum(reshape(cell2mat(BR_const_ijs_vec(:)),d*d,[]), 2); %sum all
BR_const = reshape(BR_const_vec, d, d);



end %file function

