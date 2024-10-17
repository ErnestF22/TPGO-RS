% testNetwork_test(1); %argument of som.testNetwork_test should be the test number
% run ("../manopt/importmanopt.m");
% run ("../codemeta/pathdefrepository.m");

% 1a) PW TRANSLATION DATA INPUT: R, T are the gt, Tijs_nois are the input data
% testdata = testNetwork_som(3); %4 is also the default

% som = ShapeOfMotion('testNetwork_params.csv'); %params reading is done directly in constructor
% %copy the list below from the properties list
N = 5;
d = 3;
d_aff = d+1;
global_camera_id = 1;
num_tests_per_sigma = 1;
transf_end_thresh = 1;
max_icp_iterations = 1;
num_edges_full = N*N;
num_edges = N*(N-1);
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

testdata = testNetwork_som(3);
edges = testdata.E;

Tijs_vec = G2T(testdata.gijtruth);%make_tijs_vec(points, edges, num_edges, d, N);
Tijs_mat = tijs_vec_2_tijs_mat(Tijs_vec, edges, N);
%T DATA INPUT (w.r.t. chosen camera id frame)
T_globalframe = G2T(testdata.gitruth);
R_globalframe = G2R(testdata.gitruth);

% OBS. M, N matrices are used/made inside Procrustes step 1

% 3) run Procrustes
%set initguess to truth for all rotations, and then add some random error
transf_gt = testdata.gitruth;
% R_initguess = G2R(rot_randn(testdata.gitruth, 0.0, N)); % this does not add any noise
R_initguess = R_globalframe;
transf_initguess = transf_gt;

%start timer
procrustes_start_time = tic();
rot_procrustes = som_stepone_procrustes(T_globalframe, Tijs_vec, edges, som_params);
% r_reflections = rot_procrustes(:,:,1);
% for ii = 1:N
%     rot_procrustes(:,:,ii) = rot_procrustes(:,:,ii) * r_reflections;
% end
transl_procrustes = som_steptwo_procrustes(rot_procrustes, T_globalframe, Tijs_vec, edges, som_params);
% transl_procrustes = transl_procrustes + repmat(transf_gt(1:3, 4) - transl_procrustes(1:3), N, 1);


transf_procrustes = RT2G(rot_procrustes, reshape(transl_procrustes, d, N));
% procrustes_end_time = tic();
procrustes_exec_time = toc(procrustes_start_time);

% testdata.gi = matUnstack(transf_procrustes,4);
testdata.gi = transf_procrustes;
% [rotErr,translErr,scale_ratio,translErrNorm] = testNetworkComputeErrors(testdata)
[rotErr,translErr] = testNetworkComputeErrors(testdata)

%4) compare output results

% rotation_error_procrustes = compare_rotation_results(transf_gt, procrustes_results_pairwise, som_params);
% translation_error_procrustes = compare_translation_results(transf_gt, procrustes_results_pairwise, som_params);
[rotation_error_procrustes, translation_error_procrustes] = compare_transf_results_with_gt(transf_gt, transf_procrustes, som_params);


% 5) plots of mean error and exec times

% 5a)
% plot on x sigmas,
% on y the mean (across all iterations) ||gt - result|| in terms on percentage

% 5b)
% plot on x sigmas,
% on y the mean execution time (across all iterations)
fprintf("procrustes exec time %g [s]\n", (procrustes_exec_time));

% plot (sigma, procrustes_exec_time, '.', 'markersize', 15)


