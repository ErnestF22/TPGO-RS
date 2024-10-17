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

%set data (no noise)
Tijs_vec = G2T(testdata.gijtruth);
T_globalframe = G2T(testdata.gitruth);

R_truth=G2R(testdata.gitruth);
vR_noise=rot_randTangentNormVector(R_truth);
R_initguess=rot_exp(R_truth,sigma*pi/5*vR_noise);
% transl_initguess = T_globalframe + sigma.*randn(size(T_globalframe));
% transf_initguess = RT2G(R_initguess, transl_initguess);

%add noise to data
Tijs_vec_nois = Tijs_vec + sigma.*randn(size(Tijs_vec));
T_globalframe_nois = T_globalframe + sigma.*randn(size(T_globalframe));


%% Form cost with the norm formulation
cost_norm = 0.0;
for ee = 1:size(edges,1)
    id_i = edges(ee, 1);
    id_j = edges(ee, 2);
    cost_part = norm(Tijs_vec_nois(:, ee) - ...
        R_initguess(:,:,id_i)' * T_globalframe_nois(:, id_j) + ...
        R_initguess(:,:,id_i)' * T_globalframe_nois(:, id_i))^2;
    cost_norm = cost_norm + cost_part;
end
disp("cost_norm");
disp(cost_norm);



%% Cost step1
[L, P] = make_LT_PT_noloops( ...
    T_globalframe_nois, Tijs_vec_nois, edges, som_params);
cost_const_term_tij = compute_fixed_cost_term(Tijs_vec_nois, d);

cost_step1 = trace(matStack(R_initguess)' * L * matStack(R_initguess) ...
    + matStack(R_initguess)' * P) + cost_const_term_tij;

disp("cost_step1")
disp(cost_step1)


%% Cost step2
[A, b] = make_A_b_noloops( ...
    R_initguess, T_globalframe_nois, Tijs_vec_nois, edges, som_params);

cost_step2 = norm(A*T_globalframe_nois(:) + b)^2;

disp("cost_step2")
disp(cost_step2)

