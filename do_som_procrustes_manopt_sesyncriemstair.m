function [rotation_error_manopt,translation_error_manopt, ...
    rotation_error_procrustes,translation_error_procrustes, ...
    rotation_error_sesync_riemstair,translation_error_sesync_riemstair, ...
    exectime_manopt,exectime_procrustes,exectime_sesync_riemstair] = ...
        do_som_procrustes_manopt_sesyncriemstair(testdata, sigma, mu, som_params)
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

%parse used SoM params
% N = params.N;

%edges
edges = (testdata.E);

%set gt 
% transf_gt = testdata.gitruth;

%set data (no noise)
Tijs_vec = G2T(testdata.gijtruth);
T_globalframe = G2T(testdata.gitruth);

%add noise to data
Tijs_vec_nois = Tijs_vec + sigma.*randn(size(Tijs_vec));
T_globalframe_nois = T_globalframe + sigma.*randn(size(T_globalframe));

% 3) run Manopt and then Procrustes

% 3a) execute with step 1 through MANOPT
%setup initguess
% R_initguess = G2R(rot_randn(testdata.gitruth, 0.0, N)); % this does not add any noise
R_truth=G2R(testdata.gitruth);
vR_noise=rot_randTangentNormVector(R_truth);
%R_initguess = G2R(rot_randn(testdata.gitruth, sigma_init, N));
R_initguess=rot_exp(R_truth,sigma*pi/5*vR_noise);
transl_initguess = T_globalframe + sigma.*randn(size(T_globalframe));
transf_initguess = RT2G(R_initguess, transl_initguess);
manopt_start_time = tic();
transf_manopt = som_manopt(T_globalframe_nois, Tijs_vec_nois, edges, som_params, matStack(transf_initguess));
exectime_manopt = toc(manopt_start_time);

% 3b) execute with step 1 through PROCRUSTES
procrustes_start_time = tic();
transf_procrustes = som_procrustes(T_globalframe_nois, Tijs_vec_nois, edges, som_params);
exectime_procrustes = toc(procrustes_start_time);

% 3c) execute with step 1 through Manopt with Riemannian Staircase
measurements = struct;
measurements.edges = testdata.E;
measurements.R = squeeze(mat2cell(G2R(testdata.gij), som_params.d, som_params.d, ones(som_params.num_edges, 1)))';
measurements.t = mat2cell(G2T(testdata.gij), som_params.d, ones(som_params.num_edges, 1));
measurements.kappa = num2cell(ones(1, som_params.num_edges));
measurements.tau = num2cell(ones(1, som_params.num_edges));
%SoM/testdata fields
measurements.T_globalframe_stiefel = T_globalframe_nois;
measurements.Tijs_vec = Tijs_vec_nois;
sesync_riemstair_start_time = tic();
[rot_sesync_riemstair, rot_sesync_riemstair_stief, ...
    transl_sesync_riemstair, transl_sesync_riemstair_stiefel, ...
    final_cost_sesync_riemstair, last_num_rows_stiefel] = ...
        som_riemstair_se_sync( ...
            measurements, R_initguess, transl_initguess, som_params);
exectime_sesync_riemstair = toc(sesync_riemstair_start_time);

%4) compare output results

testdata.gi = transf_manopt;
% [rotErr,translErr,scale_ratio,translErrNorm] = testNetworkComputeErrors(testdata)
[rotation_error_manopt,translation_error_manopt] = testNetworkComputeErrors(testdata);

testdata.gi = matUnstack(transf_procrustes, 4);
[rotation_error_procrustes,translation_error_procrustes] = testNetworkComputeErrors(testdata);

% testdata.gi = transf_sesync_riemstair;
testdata.gi = RT2G(rot_sesync_riemstair, transl_sesync_riemstair);
[rotation_error_sesync_riemstair,translation_error_sesync_riemstair] = testNetworkComputeErrors(testdata);


end %function

