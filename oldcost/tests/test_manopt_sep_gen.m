

% 1a) PW TRANSLATION DATA INPUT: R, T are the gt, Tijs_nois are the input data
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
riem_grad_mode = 'auto'; %'auto' or 'manual'
hessian_mode = 'auto'; %'auto' or 'manual'
initguess_is_available = boolean(1);
som_params = struct('N', N, 'd', d, 'd_aff', d_aff, ...
    'global_camera_id', global_camera_id, ...
    'num_tests_per_sigma', num_tests_per_sigma, 'transf_end_thresh', transf_end_thresh, ...
    'max_icp_iterations', max_icp_iterations, 'num_edges_full', num_edges_full, ...
    'num_edges', num_edges, 'procrustes_mode', procrustes_mode, ...
    'riem_grad_mode', riem_grad_mode, ...
    'hessian_mode', hessian_mode, ...
    'initguess_is_available', initguess_is_available);


sigma = 1.0; %change this when wanting to add noise

edges = (testdata.E);

Tijs_vec = G2T(testdata.gijtruth);
T_globalframe = G2T(testdata.gitruth);

R_initguess = G2R(rot_randn(testdata.gitruth, sigma, N));
transl_initguess = T_globalframe + sigma.*randn(size(T_globalframe));
transf_initguess.R = R_initguess;
transf_initguess.A = transl_initguess;
transf_gt = testdata.gitruth;


transf_manopt_sep = som_manopt(T_globalframe, Tijs_vec, edges, som_params, matStack(RT2G(R_initguess, transl_initguess)));

% transf_manopt_gen = som_manopt_genproc(T_globalframe, Tijs_vec, edges, som_params, matStack(transf_initguess));
transf_manopt_gen = som_manopt_genproc(T_globalframe, Tijs_vec, edges, som_params, matStack(transf_initguess));

% disp([transf_manopt_sep, transf_manopt_gen]);

%4) compare output results

testdata.gi = transf_manopt_sep;
% [rotErr,translErr,scale_ratio,translErrNorm] = testNetworkComputeErrors(testdata)
[rotation_error_manopt_sep,translation_error_manopt_sep] = testNetworkComputeErrors(testdata);

testdata.gi = transf_manopt_gen;
[rotation_error_manopt_gen,translation_error_manopt_gen] = testNetworkComputeErrors(testdata);

