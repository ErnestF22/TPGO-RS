function [rotation_error_manopt,translation_error_manopt, ...
    rotation_error_procrustes,translation_error_procrustes, ...
    rotation_error_ssom,translation_error_ssom, ...
    exectime_manopt,exectime_procrustes,exectime_ssom, ...
    scale_ratios_ssom,transl_err_norm_ssom, ...
    rs_success_bool] = ...
        do_ssom(testdata, sigma, mu, params)
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
num_edges = size(edges, 1);
testdata.edges = edges; % 2 notation for edges struct member

%% 1) add noise to data
%set gt 
% transf_gt = testdata.gitruth;

%set data (no noise)
tijs = G2T(testdata.gij);
testdata.R_gt = G2R(testdata.gitruth);
testdata.T_gt = G2T(testdata.gitruth);
testdata.lambda_gt = testdata.lambdaijtruth;
X_gt.R = testdata.R_gt;
X_gt.T = testdata.T_gt;
X_gt.lambda = testdata.lambda_gt;
testdata.tijs = tijs;
cost_gt = ssom_cost(X_gt, testdata);
disp("cost_gt in do_ssom.m")
disp(cost_gt)
% problem_data_gt.tijs = tijs;
% problem_data_gt.d = d;
% problem_data_gt.N = N;
% problem_data_gt.edges = edges;

% save('poc2degree_data/R_gt.mat', "R_globalframe")
% save('poc2degree_data/T_gt.mat', "T_globalframe")
% save('poc2degree_data/problem_data_gt.mat', "problem_data_gt")


sigma_transl = sigma;
tijs_nois = tijs + sigma_transl.*randn(size(tijs)) + ...
    mu * ones(size(tijs));

T_globalframe = G2T(testdata.gitruth);
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
T_initguess = T_globalframe + sigma.*randn(size(T_globalframe));
lambdas_initguess = ones(num_edges, 1);
if params.rand_initguess
    %overwrite sigma-noisy initguess
    R_initguess = randrot_som(params.d, params.N);
    T_initguess = 10 * rand(params.d, params.N);
    % lambdas_initguess = ones(num_edges, 1);
    % T_globalframe_nois = 10 * rand(params.d, params.N);
end
transf_initguess = RT2G(R_initguess, T_initguess);

%% 3) Run methods

% 3a) execute with step 1 through MANOPT
manopt_start_time = tic();
if params.enable_manopt_icp
    transf_manopt = rsom_manopt(T_globalframe_nois, tijs_nois, edges, params, transf_initguess);
    % manopt_end_time = tic();
else
    transf_manopt = repmat(eye(d+1), 1, 1, N);
end
exectime_manopt = toc(manopt_start_time);


% 3b) execute with step 1 through PROCRUSTES
procrustes_start_time = tic();
if params.enable_procrustes
    transf_procrustes = som_procrustes(T_globalframe_nois, tijs_nois, edges, params);
else
    transf_procrustes = repmat(eye(d+1), 1, 1, N);
end
exectime_procrustes = toc(procrustes_start_time);

% 3c) execute with step 1 through Manopt with Riemannian Staircase
ssom_start_time = tic();
if params.enable_ssom
    testdata.R_gt = X_gt.R;
    testdata.T_gt = X_gt.T;
    testdata.lambda_gt = X_gt.lambda;
    testdata.sz = [d d N];
    testdata.edges = testdata.E; %edges field name is used in rsom/ssom project, E in testnetwork benchmark testdata generator
    testdata.tijs = G2T(testdata.gij);
    [transf_ssom, lambdas_ssom_out, rs_success_bool, cost_ssom] = ...
        ssom_genproc(testdata, transf_initguess, lambdas_initguess); %lambdas_ssom_out should be used somewhere (maybe already inside ssom_genproc)
    disp("cost_ssom")
    disp(cost_ssom)
else
    rs_success_bool = boolean(0);
    transf_ssom = repmat(eye(d+1), 1, 1, N);
end
exectime_ssom = toc(ssom_start_time);

%% 4) Compare output results

testdata.gi = transf_manopt;
% [rotErr,translErr,scale_ratio,translErrNorm] = testNetworkComputeErrors(testdata)
[rotation_error_manopt,translation_error_manopt] = testNetworkComputeErrors(testdata);

testdata.gi = transf_procrustes;
[rotation_error_procrustes,translation_error_procrustes] = testNetworkComputeErrors(testdata);

testdata.gi = transf_ssom;
testdata.lambdaij = lambdas_ssom_out;
%TODO: change this back to what it should be after correcting PIM, 
%Stiefel -> SO(d) conversion
[rotation_error_ssom,translation_error_ssom,...
    scale_ratios_ssom,transl_err_norm_ssom] = ...
    testNetworkComputeErrors(testdata);
scale_ratios_ssom = (lambdas_ssom_out ./ transp(testdata.lambdaijtruth));

testNetworkDisplay(testdata);
hold on;
testNetworkDisplay(testdata, 'Estimated', 'member', 'gi');
hold off;

%rtron code
testdata=rmfield(testdata,'X')
testNetworkDisplay(testdata);
hold on;
testNetworkDisplay(testdata,'member','gi','scale',5)
hold off;

end %function

