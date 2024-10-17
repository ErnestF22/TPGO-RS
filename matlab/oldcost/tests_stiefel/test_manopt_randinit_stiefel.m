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
riem_grad_mode = 'manual'; %'auto' or 'manual'
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

transf_gt = testdata.gitruth;

som_params_subset = som_params;

progressive_transfs = zeros(d_aff,d_aff,N,N);
progressive_initguesses = zeros(d_aff,d_aff,N,N);

next_step_initguess_full.R = repmat(eye(d), 1, 1, N);
next_step_initguess_full.A = sigma.*randn(size(T_globalframe));

num_cameras_subset = 3;
subset_cameras_ids = (1:num_cameras_subset)';

first_step = boolean(1);
while num_cameras_subset <= N 
    som_params_subset.N = num_cameras_subset;
    [T_gf_subset, Tijs_vec_subset, edges_subset, som_params_subset] = ...
        make_input_subset(subset_cameras_ids, T_globalframe, Tijs_vec, edges, som_params);

    if first_step
        first_step_initguess.R = next_step_initguess_full.R(:,:,subset_cameras_ids);
        first_step_initguess.A = next_step_initguess_full.A(:,subset_cameras_ids);
        transf_initguess = first_step_initguess;
        first_step = boolean(0);
    else
        transf_initguess.R = next_step_initguess_full.R(:,:,subset_cameras_ids);
        transf_initguess.A = next_step_initguess_full.A(:,subset_cameras_ids);
        R_initguess = transf_initguess.R;
        transl_initguess = transf_initguess.A;
    end   
    progressive_initguesses(:,:,:, num_cameras_subset) = ...
        RT2G(next_step_initguess_full.R, next_step_initguess_full.A);
    
    num_rows_stiefel = d;

    transf_manopt_sep = som_manopt_stiefel(T_gf_subset, Tijs_vec_subset, edges_subset, num_rows_stiefel, som_params_subset, matStack(RT2G(transf_initguess.R, transf_initguess.A)));

    som_params.initguess_is_available = boolean(1); %for all steps after the first

    next_step_initguess_full.R(:,:,subset_cameras_ids) = G2R(transf_manopt_sep);
    next_step_initguess_full.A(:,subset_cameras_ids) = G2T(transf_manopt_sep);

    progressive_transfs(:,:,1:size(transf_manopt_sep, 3), num_cameras_subset) = ...
        transf_manopt_sep;

    num_cameras_subset = num_cameras_subset + 1;
    subset_cameras_ids = (1:num_cameras_subset)';
end

% disp([transf_manopt_sep, transf_manopt_gen]);

%4) compare output results

testdata.gi = transf_manopt_sep;
% [rotErr,translErr,scale_ratio,translErrNorm] = testNetworkComputeErrors(testdata)
[rotation_error_manopt_sep,translation_error_manopt_sep] = testNetworkComputeErrors(testdata);

% testdata.gi = transf_manopt_gen;
% [rotation_error_manopt_gen,translation_error_manopt_gen] = testNetworkComputeErrors(testdata);

