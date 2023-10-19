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

T_globalframe = G2T(testdata.gitruth);

T_globalframe_stiefel = cat_zero_row(T_globalframe);

[L_T, P_T] = make_LT_PT(T_globalframe, Tijs_vec, edges, som_params);
[L_T_stiefel, P_T_stiefel] = make_LT_PT_stiefel(T_globalframe_stiefel, Tijs_vec, edges, d_stiefel, som_params);
[L_T_nl_stiefel, P_T_nl_stiefel] = make_LT_PT_noloops_stiefel(T_globalframe_stiefel, Tijs_vec, edges, d_stiefel, som_params);

L_T = remove_quasi_zeros(L_T);
P_T = remove_quasi_zeros(P_T);
L_T_stiefel = remove_quasi_zeros(L_T_stiefel);
P_T_stiefel = remove_quasi_zeros(P_T_stiefel);
L_T_nl_stiefel = remove_quasi_zeros(L_T_nl_stiefel);
P_T_nl_stiefel = remove_quasi_zeros(P_T_nl_stiefel);

disp(max(max(abs(L_T_stiefel-L_T_nl_stiefel))))
disp(max(max(abs(P_T_stiefel-P_T_nl_stiefel))))


