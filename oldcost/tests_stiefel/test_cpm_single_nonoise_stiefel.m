% testNetwork_test(1); %argument of som.testNetwork_test should be the test number
% run ("../manopt/importmanopt.m");
% run ("../codemeta/pathdefrepository.m");

% 1a) PW TRANSLATION DATA INPUT: R, T are the gt, Tijs_nois are the input data
testdata = testNetwork_som(3); %4 is also the default

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
riem_grad_mode = 'manual'; %'auto' or 'manual'
hessian_mode = 'auto'; %'auto' or 'manual'
initguess_is_available = boolean(0);
som_params = struct('N', N, 'd', d, 'd_aff', d_aff, ...
    'global_camera_id', global_camera_id, ...
    'num_tests_per_sigma', num_tests_per_sigma, 'transf_end_thresh', transf_end_thresh, ...
    'max_icp_iterations', max_icp_iterations, 'num_edges_full', num_edges_full, ...
    'num_edges', num_edges, 'procrustes_mode', procrustes_mode, ...
    'riem_grad_mode', riem_grad_mode, ...
    'hessian_mode', hessian_mode, ...
    'initguess_is_available', initguess_is_available);


%%%

edges = (testdata.E);

Tijs_vec = G2T(testdata.gijtruth);
% Tijs_mat = tijs_vec_2_tijs_mat(Tijs_vec, edges, N);
T_globalframe = G2T(testdata.gitruth);

% OBS. M, N matrices are used/made inside Procrustes step 1

% 3) run Manopt and then Procrustes
%set initguess to truth for all rotations, and then add some random error
transf_gt = testdata.gitruth;
% R_initguess = G2R(rot_randn(testdata.gitruth, 0.0, N)); % this does not add any noise
R_initguess = G2R(testdata.gitruth);
transf_initguess = make_transf(R_initguess, T_globalframe(:));

% 3a) execute with step 1 through MANOPT
manopt_start_time = tic();
num_rows_stiefel = d;
transf_manopt = som_manopt_stiefel(T_globalframe, Tijs_vec, edges, num_rows_stiefel, som_params, transf_initguess);
% manopt_end_time = tic();
manopt_exec_time = toc(manopt_start_time);

% 3b) execute with step 1 through PROCRUSTES
procrustes_start_time = tic();
transf_procrustes = som_procrustes(T_globalframe, Tijs_vec, edges, som_params);
% procrustes_end_time = tic();
procrustes_exec_time = toc(procrustes_start_time);

%4) compare output results

testdata.gi = transf_manopt;
% [rotErr,translErr,scale_ratio,translErrNorm] = testNetworkComputeErrors(testdata)
[rotation_error_manopt,translation_error_manopt] = testNetworkComputeErrors(testdata)


testdata.gi = matUnstack(transf_procrustes, 4);
[rotation_error_procrustes,translation_error_procrustes] = testNetworkComputeErrors(testdata)

% 5) plots of mean error and exec times

% 5a)
% plot on x sigmas, 
% on y the mean (across all iterations) ||gt - result|| in terms on percentage

% 5b)
% plot on x sigmas, 
% on y the mean execution time (across all iterations) 

fprintf("manopt exec time %g [s]\n", (manopt_exec_time));

fprintf("procrustes exec time %g [s]\n", (procrustes_exec_time));

% plot (sigma, manopt_exec_time, '.', 'markersize', 15)
% plot (sigma, procrustes_exec_time, '.', 'markersize', 15)


