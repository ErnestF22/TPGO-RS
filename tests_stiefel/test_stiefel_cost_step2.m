clear;
clc;
close all;

testdata = testNetwork_som(3); %4 is the default option

% som = ShapeOfMotion('testNetwork_params.csv'); %params reading is done directly in constructor
% %copy the list below from the properties list
N = testdata.NNodes;
d = 3;
d_aff = d+1;
global_camera_id = 1;
num_tests_per_sigma = 1;
transf_end_thresh = 1;
max_icp_iterations = 1;
num_edges_full = N*N;
num_edges = testdata.NEdges;
procrustes_mode = 'som';
riem_grad_mode = 'auto'; %'auto' or 'manual'
initguess_is_available = boolean(0);
som_params = struct('N', N, 'd', d, 'd_aff', d_aff, ...
    'global_camera_id', global_camera_id, ...
    'num_tests_per_sigma', num_tests_per_sigma, 'transf_end_thresh', transf_end_thresh, ...
    'max_icp_iterations', max_icp_iterations, 'num_edges_full', num_edges_full, ...
    'num_edges', num_edges, 'procrustes_mode', procrustes_mode, ...
    'riem_grad_mode', riem_grad_mode, ...
    'initguess_is_available', initguess_is_available);


edges = (testdata.E);

d_stiefel = d+1; %d+1 = 4

Tijs_vec = G2T(testdata.gijtruth);
Tijs_mat = tijs_vec_2_tijs_mat(Tijs_vec, edges, N);

T_globalframe = G2T(testdata.gitruth);

% R_truth = G2R(testdata.gitruth);
R_truth = randrot(d, N);
R_truth_stiefel = zeros(d_stiefel, d, N);
for ii = 1:size(R_truth, 3)
    R_truth_stiefel(:,:,ii) = cat_zero_row(R_truth(:,:,ii));
end
R_truth_2d_stiefel = matStack(R_truth_stiefel);

T_globalframe_stiefel = cat_zero_row(T_globalframe);

[A, b] = make_A_b(R_truth, T_globalframe, Tijs_vec, edges, som_params);
[A_stiefel, b_stiefel] = make_A_b_stiefel(R_truth_stiefel, T_globalframe_stiefel, Tijs_vec, edges, d_stiefel, som_params);

A = remove_quasi_zeros(A);
b = remove_quasi_zeros(b);
A_stiefel = remove_quasi_zeros(A_stiefel);
b_stiefel = remove_quasi_zeros(b_stiefel);

cost_orig = A * T_globalframe(:) + b %real cost is trace(...)
cost_stiefel = A_stiefel * T_globalframe_stiefel(:) + b_stiefel %real cost is trace(...)

disp("max(cost_orig - cost_stiefel):");
disp(max(cost_orig - cost_stiefel));

