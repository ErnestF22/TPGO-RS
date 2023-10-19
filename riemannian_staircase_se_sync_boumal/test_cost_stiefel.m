close all;
clear;
clc;

%test cost in SO(d)^N

%% Create and setup test data
testdata = testNetwork_som(3); %4 would be the default

% som = ShapeOfMotion('testNetwork_params.csv'); %params reading is done directly in constructor
% %copy the list below from the properties list
N = testdata.NNodes;
d = 3;
d_aff = d+1;
global_camera_id = 1;
num_tests_per_sigma = 30;
transf_end_thresh = 1;
max_icp_iterations = 10;
num_edges_full = N*N;
num_edges = testdata.NEdges;
procrustes_mode = 'som';
riem_grad_mode = 'manual'; %'auto' or 'manual'
hessian_mode = 'auto'; 
initguess_is_available = boolean(1);
som_params = struct('N', N, 'd', d, 'd_aff', d_aff, ...
    'global_camera_id', global_camera_id, ...
    'num_tests_per_sigma', num_tests_per_sigma, 'transf_end_thresh', transf_end_thresh, ...
    'max_icp_iterations', max_icp_iterations, 'num_edges_full', num_edges_full, ...
    'num_edges', num_edges, 'procrustes_mode', procrustes_mode, ...
    'riem_grad_mode', riem_grad_mode, ...
    'hessian_mode', hessian_mode, ...
    'initguess_is_available', initguess_is_available);

sigma = 0.0;

edges = (testdata.E);

num_rows_stiefel = d+1;

%set data (no noise)
Tijs_vec = G2T(testdata.gijtruth);
T_globalframe = G2T(testdata.gitruth);

R_truth=G2R(testdata.gitruth);
vR_noise=rot_randTangentNormVector(R_truth);
R_initguess=rot_exp(R_truth,sigma*pi/5*vR_noise);
R_initguess_stiefel = cat_zero_cols_3d_array(R_initguess, num_rows_stiefel-d);

%add noise to data
Tijs_vec_nois = Tijs_vec + sigma.*randn(size(Tijs_vec));
T_globalframe_nois = T_globalframe + sigma.*randn(size(T_globalframe));
T_globalframe_nois_stiefel = cat_zero_row(T_globalframe_nois, num_rows_stiefel-d);

Y_initguess_stiefel = multitransp(R_initguess_stiefel);
transl_initguess = T_globalframe_nois;
% transf_initguess = RT2G(multitransp(R_initguess), transl_initguess);

%% Form cost with the norm formulation
cost_norm = 0.0;
for ee = 1:size(edges,1)
    id_i = edges(ee, 1);
    id_j = edges(ee, 2);
    cost_part = norm(Tijs_vec_nois(:, ee) - ...
        Y_initguess_stiefel(:,:,id_i)' * T_globalframe_nois_stiefel(:, id_j) + ...
        Y_initguess_stiefel(:,:,id_i)' * T_globalframe_nois_stiefel(:, id_i))^2;
    cost_norm = cost_norm + cost_part;
end
disp("cost_norm");
disp(cost_norm);

cost_const_term_tij = compute_fixed_cost_term(Tijs_vec_nois, d);

transf_out = som_stepone_manopt_stiefel_boumal( ...
    T_globalframe_nois, Tijs_vec_nois, edges, cost_const_term_tij, ...
    num_rows_stiefel, som_params, Y_initguess_stiefel);



%% Cost step1
[L_stiefel, P_stiefel] = make_LT_PT_noloops_stiefel( ...
    T_globalframe_nois_stiefel, Tijs_vec_nois, edges, num_rows_stiefel, som_params);
cost_const_term_tij = compute_fixed_cost_term(Tijs_vec_nois, d);

cost_step1 = trace(matStack(Y_initguess_stiefel)' * L_stiefel * matStack(Y_initguess_stiefel) ...
    + matStack(Y_initguess_stiefel)' * P_stiefel) + cost_const_term_tij;

disp("cost_step1")
disp(cost_step1)


%% Cost step2
[A_stiefel, b_stiefel] = make_A_b_noloops_stiefel( ...
    Y_initguess_stiefel, T_globalframe_nois_stiefel, Tijs_vec_nois, edges, ...
    num_rows_stiefel, som_params);

cost_step2 = norm(A_stiefel*T_globalframe_nois_stiefel(:) + b_stiefel)^2;

disp("cost_step2")
disp(cost_step2)

