% testNetwork_test(1); %argument of som.testNetwork_test should be the test number
% run ("../manopt/importmanopt.m");
% run ("../codemeta/pathdefrepository.m");

% 0a) SOM PARAMS
som = ShapeOfMotion('noise_test_params.csv'); %params reading is done directly in constructor
%copy the list below from the properties list
N = som.N;
d = som.d;
d_aff = som.d_aff;
global_camera_id = som.global_camera_id;
num_tests_per_sigma = som.num_tests_per_sigma;
transf_end_thresh = som.transf_end_thresh;
max_icp_iterations = som.max_icp_iterations;
num_edges_full = som.num_edges_full;
num_edges = som.num_edges;
procrustes_mode = som.procrustes_mode;
initguess_is_available = som.initguess_is_available;
% initguess_is_available = boolean(0);
som_params = struct('N', N, 'd', d, 'd_aff', d_aff, ...
    'global_camera_id', global_camera_id, ...
    'num_tests_per_sigma', num_tests_per_sigma, 'transf_end_thresh', transf_end_thresh, ...
    'max_icp_iterations', max_icp_iterations, 'num_edges_full', num_edges_full, ...
    'num_edges', num_edges, 'procrustes_mode', procrustes_mode, ...
    'initguess_is_available', initguess_is_available);

% 0b) SOM PARAMS
%NOTE: sigmas, mus can be seen as couples for each test
sigmas = readmatrix("sigmas.txt"); %sigma = stdev, sigma.^2 = variance
mus = readmatrix("../data/mus.txt"); %OBS. generally, mus can be d-dimensional; here, we just assume them as scalar (i.e. a d-dimensional vector with all coordinates equal)

%TOUSE: when multple sigmas and mus, and multiple tests per each pair
% for ii = 1:size(sigmas,1)
%     noise_params = struct('sigma', sigmas(ii), 'mu', mus(ii)); 
%     do_som();
% end    
noise_params = struct('sigma', sigmas(1), 'mu', mus(1)); 
sigma = noise_params.sigma;
mu = noise_params.mu;

%%%

edges = (testdata.E);

Tijs_vec = G2T(testdata.gijtruth);
Tijs_mat = tijs_vec_2_tijs_mat(Tijs_vec, edges, d, N);

T_globalframe = G2T(testdata.gitruth);

T_globalframe_nois = setup_noise_tests(T_globalframe, sigma, mu);
% R_globalframe_nois = setup_noise_tests(R_globalframe, sigma, mu);

% OBS. M, N matrices are used/made inside Procrustes step 1

% 3) run Manopt and then Procrustes
%set initguess to truth for all rotations
R_initguess = G2R(testdata.gitruth);
transf_initguess = make_transf(R_initguess, T_globalframe(:));

% 3a) execute with step 1 through MANOPT
manopt_start_time = tic();
transf_manopt = som_simple_coord_desc_manopt(T_globalframe_nois, Tijs_vec, Tijs_mat, edges, L_T, P_T, som_params, transf_initguess);
% manopt_end_time = tic();
manopt_exec_time = toc(manopt_start_time);

% 3b) execute with step 1 through PROCRUSTES
procrustes_start_time = tic();
transf_procrustes = som_simple_coord_desc_procrustes(T_globalframe_nois, Tijs_vec, Tijs_mat, edges, som_params);
% procrustes_end_time = tic();
procrustes_exec_time = toc(procrustes_start_time);


% for converting affine transf matrices into rigid transform objects, use:
% rigidtform3d(transf_procrustes(1:4,:))


% 4) plots

% 4a)
% plot on x sigmas, 
% on y the mean (across all iterations) ||gt - result|| in terms on percentage

%make transf_gt 
transf_gt = R_initguess; %first, make it with rotation only, then expand it
for ii = 1:N
    transf_gt(d_aff,d_aff,ii) = 0; %first, expand with zeros
    transf_gt(d_aff,d_aff,ii) = 1; %set bottom right (affine scale) element to 1
    transf_gt(1:d,d_aff,ii) = T_globalframe(:, ii); %set equal to Tijs (gt)
end

procrustes_results_pairwise = make_results_pairwise(transf_procrustes, som_params);
manopt_results_pairwise = make_results_pairwise(transf_manopt, som_params); 

% rotation_error_manopt = compare_rotation_results(transf_gt, manopt_results_pairwise, som_params);
% translation_error_manopt = compare_translation_results(transf_gt, manopt_results_pairwise, som_params);
[rotation_error_manopt, translation_error_manopt] = compare_transf_results(transf_gt, manopt_results_pairwise, som_params);

% rotation_error_procrustes = compare_rotation_results(transf_gt, procrustes_results_pairwise, som_params);
% translation_error_procrustes = compare_translation_results(transf_gt, procrustes_results_pairwise, som_params);
[rotation_error_procrustes, translation_error_procrustes] = compare_transf_results(transf_gt, procrustes_results_pairwise, som_params);


% 4b)
% plot on x sigmas, 
% on y the mean execution time (across all iterations) 

disp("manopt_exec_time");
disp(manopt_exec_time);

disp("procrustes_exec_time");
disp(procrustes_exec_time);

manopt_exec_times = manopt_exec_time .* ones(size(sigmas)); %FIXME: substitute ones(size(sigmas)) with sigmas when tests with multiple sigmas implemented
plot (sigmas, manopt_exec_times, '.', 'markersize', 15)



