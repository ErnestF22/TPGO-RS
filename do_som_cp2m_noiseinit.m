function [rotation_error_som_riemstair,translation_error_som_riemstair, ...
    rotation_error_procrustes,translation_error_procrustes, ...
    rotation_error_manopt_sep,translation_error_manopt_sep, ...
    rotation_error_manopt_gen,translation_error_manopt_gen, ...
    exectime_som_riemstair, ...
    exectime_procrustes, ...
    exectime_manopt, ...
    exectime_manopt_sep] = do_som_cp2m_noiseinit(testdata, ...
    sigma_noise, sigma_init, mu, params)
%DO_SOM_CP2M_NOISEINIT Function that executes the Shape of Motion 
%algorithms through Procrustes, "Separated" Manopt and 
%"Generalized Procrustes" Manopt, returning  rotation and translation 
%errors, as well as the mean execution time;
%the SoM algorithms are run on Gaussian noisy input data, 
%with sigma_noise variance on the Tijs inputs and sigma_init noise on the 
%initial guesses of rotation and translation and mu mean (that are being 
%passed as arguments to the DO_SOM function)
%The generalized Procrustes pipeline for Manopt optimizes at the same time 
%for rotations and translations, opposite to the do_som() original two-steps pipeline.

%parse used SoM params
N = params.N;

%edges
edges = (testdata.E);

%set gt 
% transf_gt = testdata.gitruth;

%set data (no noise)
Tijs_vec = G2T(testdata.gijtruth);
T_globalframe = G2T(testdata.gitruth);

%add noise to data
Tijs_vec_nois = Tijs_vec + sigma_noise.*randn(size(Tijs_vec));
%T_globalframe_nois = T_globalframe + sigma_init.*randn(size(T_globalframe));
T_noise=cnormalize(randn(size(T_globalframe)));
T_globalframe_nois = T_globalframe + sigma_init*T_noise;

% 3) run Manopt and then Procrustes

% 3a) execute with step 1 through MANOPT SEP
%setup initguess
% R_initguess = G2R(rot_randn(testdata.gitruth, 0.0, N)); % this does not add any noise
R_truth=G2R(testdata.gitruth);
vR_noise=rot_randTangentNormVector(R_truth);
%R_initguess = G2R(rot_randn(testdata.gitruth, sigma_init, N));
R_initguess=rot_exp(R_truth,sigma_init*pi/5*vR_noise);
transl_initguess = T_globalframe_nois;
transf_initguess.R = R_initguess;
transf_initguess.A = transl_initguess;
% transf_gt = testdata.gitruth;

% 3
% 3a) execute with step 1 through SOM_RIEMANNIAN_STAIRCASE
som_riemstair_start_time = tic();
transf_som_riemstair = som_riemannian_staircase(T_globalframe_nois, Tijs_vec_nois, edges, params, R_initguess, transl_initguess);
% procrustes_end_time = tic();
exectime_som_riemstair = toc(som_riemstair_start_time);

% 3b) execute with step 1 through PROCRUSTES
procrustes_start_time = tic();
transf_procrustes = som_procrustes(T_globalframe_nois, Tijs_vec_nois, edges, params);
% procrustes_end_time = tic();
exectime_procrustes = toc(procrustes_start_time);

% 3c) execute with step 1 through MANOPT
manopt_start_time = tic();
transf_manopt_sep = som_manopt(T_globalframe_nois, Tijs_vec_nois, edges, params, matStack(RT2G(R_initguess, transl_initguess)));
% manopt_end_time = tic();
exectime_manopt = toc(manopt_start_time);

% 3d) execute with step 1 through MANOPT GENPROC
manopt_sep_start_time = tic();
transf_manopt_gen = som_manopt_genproc(T_globalframe_nois, Tijs_vec_nois, edges, params, matStack(transf_initguess));
% procrustes_end_time = tic();
exectime_manopt_sep = toc(manopt_sep_start_time);


%4) compare output results
% [rotErr,translErr,scale_ratio,translErrNorm] = testNetworkComputeErrors(testdata)

testdata.gi = transf_som_riemstair;
[rotation_error_som_riemstair,translation_error_som_riemstair] = ...
    testNetworkComputeErrors(testdata);

testdata.gi = matUnstack(transf_procrustes, 4);
[rotation_error_procrustes,translation_error_procrustes] = ...
    testNetworkComputeErrors(testdata);

testdata.gi = transf_manopt_sep;
[rotation_error_manopt_sep,translation_error_manopt_sep] = ...
    testNetworkComputeErrors(testdata);

testdata.gi = transf_manopt_gen;
[rotation_error_manopt_gen,translation_error_manopt_gen] = ...
    testNetworkComputeErrors(testdata);

end

