function P = make_p_noloops(R_gf, T_gf, Tijs_vec, edges)
%MAKE_P_NOLOOPS Make matrix P used in step1 (rotation estimation) 
% cost formulation WITHOUT using any for loops
nrs = size(R_gf, 1);
d = size(R_gf, 2);
N = size(R_gf, 3);

adj_mat = make_adj_mat_from_edges(edges, N);

M = num2cell(adj_mat);
% P = zeros(nrs, d*N);

idx_col_p = reshape(1:d*N, [], N)';

num_edges = size(edges,1);

%A
T_gf_cells = num2cell(T_gf, 1);
A_1 = repmat(T_gf_cells, N, 1);
A_2 = repmat(T_gf_cells', 1, N);
A = cellfun(@minus,A_1,A_2,'Un',0);

%
Tijs_mat = tijs_vec_2_tijs_mat(Tijs_vec, edges, N);
B = mat2cell(Tijs_mat, d*ones(1,N), ones(1,N));
B_transp = cellfun(@transp,B,'Un',0);

C = cellfun(@times,A,B_transp,'Un',0);
% celldisp(C)

P_half = cellfun(@times,C,M,'Un',0); %apply mask
P = 2 * cell2mat(P_half);

end %file function


