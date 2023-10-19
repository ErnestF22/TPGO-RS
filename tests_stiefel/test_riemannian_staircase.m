clear;
clc;
close all;

testdata = testNetwork_som(3); %4 would be the default


% 0a) SOM PARAMS
% som = ShapeOfMotion('noise_test_params.csv'); %params reading is done directly in constructor
%copy the list below from the properties list
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
hessian_mode = 'manual'; %'auto' or 'manual'
initguess_is_available = boolean(0);
som_params = struct('N', N, 'd', d, 'd_aff', d_aff, ...
    'global_camera_id', global_camera_id, ...
    'num_tests_per_sigma', num_tests_per_sigma, 'transf_end_thresh', transf_end_thresh, ...
    'max_icp_iterations', max_icp_iterations, 'num_edges_full', num_edges_full, ...
    'num_edges', num_edges, 'procrustes_mode', procrustes_mode, ...
    'riem_grad_mode', riem_grad_mode, ...
    'hessian_mode', hessian_mode, ...
    'initguess_is_available', initguess_is_available);


% 0b) SOM PARAMS
%NOTE: sigmas, mus can be seen as couples for each test
sigmas = readmatrix("sigmas.txt"); %sigma = stdev, sigma.^2 = variance
mus = readmatrix("mus.txt"); %OBS. generally, mus can be d-dimensional; here, we just assume them as scalar (i.e. a d-dimensional vector with all coordinates equal)

%TOUSE: when multple sigmas and mus, and multiple tests per each pair
% for ii = 1:size(sigmas,1)
%     noise_params = struct('sigma', sigmas(ii), 'mu', mus(ii)); 
%     do_som();
% end    
noise_params = struct('sigma', sigmas(1), 'mu', mus(1)); 
sigma = noise_params.sigma;
mu = noise_params.mu;


% 1a) PW TRANSLATION DATA INPUT: R, T are the gt, Tijs_mat_nois are the input data
%Rijs, Tijs are the used just for making these dummy tests without having
%the real data

edges = (testdata.E); %graph is not full by default

Tijs_vec = G2T(testdata.gijtruth);
T_globalframe = G2T(testdata.gitruth);

Tijs_mat = tijs_vec_2_tijs_mat(Tijs_vec, edges, N);
%T DATA INPUT (w.r.t. chosen camera id frame)

R_initguess = G2R(testdata.gitruth);
transl_initguess = T_globalframe;
transf_initguess.R = R_initguess;
transf_initguess.A = transl_initguess;
transf_gt = testdata.gitruth;

T_globalframe_nois = T_globalframe;
Tijs_vec_nois = Tijs_vec;
Tijs_mat_nois = Tijs_mat;

% transf_curr = RT2G(R, reshape(T, d, N));
% disp(transf_curr);

transf_out = som_riemannian_staircase( ...
    T_globalframe, Tijs_vec, edges, som_params, R_initguess, transl_initguess);

% [R_out, T_out] = lowRankLocalization_solution_extractProjection(matStack(multitransp(Y_star_riemstaircase)) * reshape(T_stiefel, [], N));

Y_star_riemstaircase = G2R(transf_out);
T_stiefel = G2T(transf_out);

disp("Y_star_riemstaircase");
disp(Y_star_riemstaircase);
disp("T_stiefel");
disp(T_stiefel);


