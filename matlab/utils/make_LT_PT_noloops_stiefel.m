function [L_T_noloops, P_T_noloops] = make_LT_PT_noloops_stiefel(T_globalframe_nois, Tijs_vec, edges, d_stiefel, params)
%MAKE_LT_PT L_T will be a block-diagonal matrix, with i-th block equal to the sum of 4
%terms (L_one - L_two - L_three + L_four)
%L_one_i = sum_j T_i*T_j'
%L_two_i = sum_i T_i*T_i'
%L_three_i = sum_j T_j*T_j'
%L_four_i = sum_j T_j*T_i'
%P_T will be dNxd matrix, with i-th block-row equal to:
% \sum_{j=1}^{N} P_{ij}
% with P_{ij} = 2(T_i-T_j)T_{ij}\transpose

N = params.N;
d = params.d;
d_aff = params.d_aff;
global_camera_id = params.global_camera_id;
num_tests_per_sigma = params.num_tests_per_sigma;
transf_end_thresh = params.transf_end_thresh;
max_icp_iterations = params.max_icp_iterations;
num_edges_full = params.num_edges_full;
num_edges = params.num_edges;
initguess_is_available = params.initguess_is_available;

% L_T_noloops = zeros(d*N, d*N); %would be overwritten anyways
% P_T_noloops = zeros(d*N,d); %would be overwritten anyways

%make adjacency, incidence, directed incidence matrices
incidence_mat = make_incidence_matrix_from_edges(edges, N);
incidence_mat_dir = make_directed_incidence_matrix_from_edges(edges, N);
adj_mat = make_adj_mat_from_edges(edges, N);
adj_mat_3d_colmaj = make_adj_mat_3d_colmaj(adj_mat);

Tijs_mat = tijs_vec_2_tijs_mat(Tijs_vec, edges, N); %add to overleaf

Tijs_3d_colmaj = make_tijs_3d(Tijs_mat); %used only for P_T

T_globalframe_3d = make_t_globalframe_3d(T_globalframe_nois);
T_globalframe_3d_repelem = repelem(T_globalframe_3d, 1,1,N);
T_globalframe_3d_repmat = repmat(T_globalframe_3d, 1,1,N);

%fill L(T)
T_globalframe_exp_orig = multiprod(T_globalframe_3d_repelem, multitransp(T_globalframe_3d_repmat)); %exp stands for expanded
shape_orig = size(T_globalframe_exp_orig);
T_globalframe_exp_orig_weighted = multiprod(adj_mat_3d_colmaj, T_globalframe_exp_orig);

%L_1
T_globalframe_exp_one = reshape(T_globalframe_exp_orig_weighted, d_stiefel, d_stiefel, N, N);
L_one = sum(T_globalframe_exp_one, 3);
% L_one_blkdiag =
% blkdiag(L_one(:,1:3),L_one(:,4:6),L_one(:,7:9),L_one(:,10:12),L_one(:,13:15))
% %works only with N=5
L_one_cells = mat2cell(matStack(squeeze(L_one)), ones(N,1)*d_stiefel, d_stiefel);
L_one_blkdiag = blkdiag(L_one_cells{:});

%L_2
% mask_two_blkdiag = kron(eye(N),ones(d));
T_two = multiprod(T_globalframe_3d, multitransp(T_globalframe_3d));
imat_dir_two = incidence_mat_dir;
imat_dir_two(imat_dir_two>=0) = 0;
imat_dir_two_colsum = - sum(imat_dir_two, 2); %!! - sum
imat_dir_two_colsum_3d = reshape(imat_dir_two_colsum, 1, 1, N);
L_two_unstacked = multiprod(imat_dir_two_colsum_3d, T_two);
L_two_cells = mat2cell(matStack(L_two_unstacked), ones(N,1)*d_stiefel, d_stiefel);
L_two_blkdiag = blkdiag(L_two_cells{:});

%L_3
T_three = multiprod(T_globalframe_3d, multitransp(T_globalframe_3d));
T_three_rep = repmat(T_three, 1, 1, N);
adj_mat_zerodiag = adj_mat;
adj_mat_zerodiag(1:1+size(adj_mat_zerodiag, 1):end) = 0; %set all diag elements to 0 (just in case)
adj_mat_zerodiag_3d_colmaj = make_adj_mat_3d_colmaj(adj_mat_zerodiag);
T_globalframe_exp_three = multiprod(adj_mat_zerodiag_3d_colmaj, T_three_rep);
L_3_to_be_summed = reshape(T_globalframe_exp_three, d_stiefel, d_stiefel, N, N);
L_3_summed = sum(L_3_to_be_summed, 3);
L_three_cells = mat2cell(matStack(squeeze(L_3_summed)), ones(N,1)*d_stiefel, d_stiefel);
L_three_blkdiag = blkdiag(L_three_cells{:});


%L_4
T_globalframe_exp_four = reshape(T_globalframe_exp_orig_weighted, d_stiefel, d_stiefel, N, N);
L_four = sum(T_globalframe_exp_four, 4);
L_four_cells = mat2cell(matStack(squeeze(L_four)), ones(N,1)*d_stiefel, d_stiefel);
L_four_blkdiag = blkdiag(L_four_cells{:});

L_T_noloops = - L_one_blkdiag + L_two_blkdiag + L_three_blkdiag - L_four_blkdiag;


%fill P(T)
P_globalframe_exp_orig = T_globalframe_3d_repmat - T_globalframe_3d_repelem;
P_T_noloops_unweighted_exp = 2.* multiprod(P_globalframe_exp_orig, multitransp(Tijs_3d_colmaj));
P_T_noloops_weighted_exp = multiprod(adj_mat_3d_colmaj, P_T_noloops_unweighted_exp); %weighted
P_T_noloops_weighted_exp_4d = reshape(P_T_noloops_weighted_exp, d_stiefel, d, N, N);
P_T_noloops_4d = sum(P_T_noloops_weighted_exp_4d, 4);
P_T_noloops = matStack(squeeze(P_T_noloops_4d));

end

