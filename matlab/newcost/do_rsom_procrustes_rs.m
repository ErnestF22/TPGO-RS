function [rotation_error_manopt,translation_error_manopt, ...
    rotation_error_procrustes,translation_error_procrustes, ...
    rotation_error_manopt_rs,translation_error_manopt_rs, ...
    exectime_manopt,exectime_procrustes,exectime_manopt_rs] = ...
        do_rsom_procrustes_rs(testdata, sigma, mu, params)
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
sigma_transl = sigma;
Tijs_vec_nois = Tijs_vec + sigma_transl.*randn(size(Tijs_vec)) + ...
    mu * ones(size(Tijs_vec));
T_globalframe_nois = T_globalframe + sigma_transl.*randn(size(T_globalframe)) + ...
    mu * ones(size(T_globalframe));

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
transf_manopt = rsom_manopt(T_globalframe_nois, Tijs_vec_nois, edges, params, transf_initguess);
% manopt_end_time = tic();
exectime_manopt = toc(manopt_start_time);

% 3b) execute with step 1 through PROCRUSTES
procrustes_start_time = tic();
transf_procrustes = som_procrustes(T_globalframe_nois, Tijs_vec_nois, edges, params);
% procrustes_end_time = tic();
exectime_procrustes = toc(procrustes_start_time);

% 3c) execute with step 1 through Manopt with Riemannian Staircase
manopt_rc_start_time = tic();
transf_manopt_rs = rsom_rs(T_globalframe_nois, Tijs_vec_nois, edges, params, transf_initguess);
% manopt_rc_end_time = tic();
exectime_manopt_rs = toc(manopt_rc_start_time);

%4) compare output results

testdata.gi = transf_manopt;
% [rotErr,translErr,scale_ratio,translErrNorm] = testNetworkComputeErrors(testdata)
[rotation_error_manopt,translation_error_manopt] = testNetworkComputeErrors(testdata);

testdata.gi = matUnstack(transf_procrustes, 4);
[rotation_error_procrustes,translation_error_procrustes] = testNetworkComputeErrors(testdata);

testdata.gi = transf_manopt_rs;
%TODO: change this back to what it should be after correcting PIM, 
%Stiefel -> SO(d) conversion
[rotation_error_manopt_rs,translation_error_manopt_rs] = testNetworkComputeErrors(testdata);


end %function

