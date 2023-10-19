function [A_nl,b_nl,R_stiefel] = make_A_b_noloops_stiefel(R_stiefel, T_gf_stiefel, Tijs_vec, edges, d_stiefel, params)
%MAKE_A_B function returns A,b matrices used in step 2 of Manopt-related
%translation estimation procedure
%and covert R to "3D" mode in case it is passed as stacked

N = params.N;
d = params.d;
d_aff = params.d_aff;
global_camera_id = params.global_camera_id;
num_tests_per_sigma = params.num_tests_per_sigma;
transf_end_thresh = params.transf_end_thresh;
max_icp_iterations = params.max_icp_iterations;
num_edges_full = params.num_edges_full;
num_edges = params.num_edges;
procrustes_mode = params.procrustes_mode;
initguess_is_available = params.initguess_is_available;



%make adjacency, incidence, directed incidence matrices
incidence_mat = make_incidence_matrix_from_edges(edges, N);
incidence_mat_dir = make_directed_incidence_matrix_from_edges(edges, N);
adj_mat = make_adj_mat_from_edges(edges, N);
adj_mat_3d_colmaj = make_adj_mat_3d_colmaj(adj_mat);


%A_nl
%A(col_id*d-(d-1):col_id*d, (ii*d)-(d-1):ii*d) = R(:,:,ii)';
%A(col_id*d-(d-1):col_id*d, (jj*d)-(d-1):jj*d) = -R(:,:,ii)';
R_multitransp = multitransp(R_stiefel);
R_multitransp_repelem = repelem(R_multitransp, 1, 1, num_edges);
R_multitransp_repmat = repmat(R_multitransp, 1, 1, num_edges);
% R_multitransp_h_stack = matStackH_d2(R_multitransp);
% R_multitransp_h_stack_rep = repmat(R_multitransp_h_stack, num_edges, 1);
incidence_mat_dir_transp = incidence_mat_dir';
incidence_mat_3d_colmaj = reshape(incidence_mat_dir_transp, 1, 1, []);
out_indices = find(incidence_mat_3d_colmaj>0);
in_indices = find(incidence_mat_3d_colmaj<0);
incidence_mat_3d_colmaj_in = incidence_mat_3d_colmaj;
incidence_mat_3d_colmaj_in(:,:,out_indices) = 0;
incidence_mat_3d_colmaj_out = incidence_mat_3d_colmaj;
incidence_mat_3d_colmaj_out(:,:,in_indices) = 0;

Anl_3d = multiprod(R_multitransp_repelem, - incidence_mat_3d_colmaj_in);
% Anl_3d(:,:,out_indices) = Anl_3d(:,:,in_indices);
Anl_4d = reshape(Anl_3d, d_stiefel, d, num_edges, N);
Anl_4d_perm = permute(Anl_4d, [1,2,4,3]);
Anl_hor_stacked = matStackH_d2(Anl_4d_perm);
Anl_hor_stacked_cells = mat2cell(Anl_hor_stacked, d, (d_stiefel*N).*ones(1, num_edges));
A_nl = cell2mat(Anl_hor_stacked_cells(:));

A_nl_copy = A_nl;
A_nl_copy_cells = mat2cell(A_nl_copy, d*ones(1, num_edges), d_stiefel*ones(1, N));
incidence_mat_dir_out = incidence_mat_dir';
incidence_mat_dir_out(incidence_mat_dir_out>0)=0;
incidence_mat_dir_out = - incidence_mat_dir_out;
incidence_mat_dir_in = incidence_mat_dir';
incidence_mat_dir_in(incidence_mat_dir_in<0)=0;
find_out = (incidence_mat_dir_out>0);
find_in = (incidence_mat_dir_in>0);
find_out_transp = find_out';
find_in_transp = find_in';
A_nl_copy_cells_transp = A_nl_copy_cells';
A_nl_copy_cells_transp(find_in_transp) = A_nl_copy_cells_transp(find_out_transp); 
A_nl_copy_cells_transp(find_in_transp) = cellfun(@(x) x*-1,A_nl_copy_cells_transp(find_in_transp),'un',0);
A_nl_copy_cells = A_nl_copy_cells_transp';
A_nl = cell2mat(A_nl_copy_cells);

% A_nl = cell2mat(reshape(A_nl_out_copy_cells_transp_vec, num_edges, N));

%b_nl
b_nl = Tijs_vec(:);

end %function
