function [rotation_error_manopt,translation_error_manopt, ...
    rotation_error_procrustes,translation_error_procrustes, ...
    rotation_error_manopt_genproc,translation_error_manopt_genproc, ...
    exectime_manopt,exectime_procrustes,exectime_manopt_genproc, ...
    rs_success_bool] = ...
        do_rsom_procrustes_genproc(testdata, sigma, mu, params)
%DO_SOM_PROCRUSTES_MANOPT_RIEMANNIAN_STAIRCASE
%Function that executes the Shape of Motion algorithms through
%Manopt and Procrustes pipelines, as well as the Manopt with the added 
%Riemannian Staircase ICP, returning rotation and translation
%errors, as well as the mean execution time; the SoM algorithms are run on
%Gaussian noisy input data, with sigma variance and mu mean (that are being 
%passed as arguments to the DO_SOM function)
%This version uses the two steps pipeline for Manopt, where optimization is
%done repetitively first on rotations and then on translations (in an
%ICP style pipeline).

if ~exist('mu', 'var')
    mu = 0.0;
end


%% 0) parse used SoM params
N = params.N;
d = params.d;

%edges
edges = (testdata.E);

%% 1) add noise to data
%set gt 
% transf_gt = testdata.gitruth;

%set data (no noise)
Tijs_vec = G2T(testdata.gijtruth);
T_globalframe = G2T(testdata.gitruth);
R_globalframe = G2R(testdata.gitruth);
X_gt.R = R_globalframe;
X_gt.T = T_globalframe;
problem_data_gt.Tijs = Tijs_vec;
problem_data_gt.d = 3;
problem_data_gt.N = 6;
problem_data_gt.edges = edges;
disp("rsom_cost_base GT")
disp(rsom_cost_base(X_gt, problem_data_gt))
% save('poc2degree_data/R_gt.mat', "R_globalframe")
% save('poc2degree_data/T_gt.mat', "T_globalframe")
% save('poc2degree_data/problem_data_gt.mat', "problem_data_gt")


sigma_transl = sigma;
Tijs_vec_nois = Tijs_vec + sigma_transl.*randn(size(Tijs_vec)) + ...
    mu * ones(size(Tijs_vec));

T_globalframe_nois = T_globalframe + sigma_transl.*randn(size(T_globalframe)) + ...
    mu * ones(size(T_globalframe));

if sigma == 0
    params.noisy_test = boolean(0);
else
    params.noisy_test = boolean(1);
end


%% 2) setup initguess
% R_initguess = G2R(rot_randn(testdata.gitruth, 0.0, N)); % this does not add any noise
R_truth=G2R(testdata.gitruth);
vR_noise=rot_randTangentNormVector(R_truth);
%R_initguess = G2R(rot_randn(testdata.gitruth, sigma_init, N));
R_initguess=rot_exp(R_truth,sigma*pi/5*vR_noise);
transl_initguess = T_globalframe + sigma.*randn(size(T_globalframe));
if params.rand_initguess
    %overwrite sigma-noisy initguess
    R_initguess = randrot_som(params.d, params.N);
    transl_initguess = 10 * rand(params.d, params.N);
    %
%     T_globalframe_nois = 10 * rand(params.d, params.N);
end
transf_initguess = RT2G(R_initguess, transl_initguess);

%% 3) Run methods

% 3a) execute with step 1 through MANOPT
manopt_start_time = tic();
if params.enable_manopt_icp
    transf_manopt = rsom_manopt(T_globalframe_nois, Tijs_vec_nois, edges, params, transf_initguess);
    % manopt_end_time = tic();
else
    transf_manopt = repmat(eye(d+1), 1, 1, N);
end
exectime_manopt = toc(manopt_start_time);


% 3b) execute with step 1 through PROCRUSTES
procrustes_start_time = tic();
if params.enable_procrustes
    transf_procrustes = som_procrustes(T_globalframe_nois, Tijs_vec_nois, edges, params);
else
    transf_procrustes = repmat(eye(d+1), 1, 1, N);
end
exectime_procrustes = toc(procrustes_start_time);

% 3c) execute with step 1 through Manopt with Riemannian Staircase
manopt_genproc_start_time = tic();
if params.enable_manopt_rs
    params.R_gt = X_gt.R;
    params.T_gt = X_gt.T;
    params.testdata = testdata; % !!
    [transf_manopt_rs, rs_success_bool, cost_manopt_rs] = ...
        rsom_genproc(T_globalframe_nois, Tijs_vec_nois, edges, params, transf_initguess);
    disp("cost_manopt_rs")
    disp(cost_manopt_rs)
else
    rs_success_bool = boolean(0);
    transf_manopt_rs = repmat(eye(d+1), 1, 1, N);
end
exectime_manopt_genproc = toc(manopt_genproc_start_time);

%% 4) Compare output results

testdata.gi = transf_manopt;
% [rotErr,translErr,scale_ratio,translErrNorm] = testNetworkComputeErrors(testdata)
[rotation_error_manopt,translation_error_manopt] = testNetworkComputeErrors(testdata);

testdata.gi = transf_procrustes;
[rotation_error_procrustes,translation_error_procrustes] = testNetworkComputeErrors(testdata);

testdata.gi = transf_manopt_rs;
%TODO: change this back to what it should be after correcting PIM, 
%Stiefel -> SO(d) conversion
[rotation_error_manopt_genproc,translation_error_manopt_genproc] = ...
    testNetworkComputeErrors(testdata);


end %function

