function [rotation_error_manopt_sep,translation_error_manopt_sep, ...
    rotation_error_manopt_gen,translation_error_manopt_gen, ...
    exectime_manopt,exectime_manopt_sep] = do_som_genproc(testdata, sigma, mu, params)
%DO_SOM_GENPROC Function that executes the Shape of Motion algorithms through
%"Separated" Manopt and "Generalized Procrustes" Manopt, returning 
%rotation and translation errors, as well as the mean execution time;
%the SoM algorithms are run on
%Gaussian noisy input data, with sigma variance and mu mean (that are being 
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
transf_initguess.R = R_initguess;
transf_initguess.A = transl_initguess;
% transf_gt = testdata.gitruth;
manopt_start_time = tic();
transf_manopt_sep = som_manopt(T_globalframe_nois, Tijs_vec_nois, edges, params, matStack(RT2G(R_initguess, transl_initguess)));
% manopt_end_time = tic();
exectime_manopt = toc(manopt_start_time);

% 3b) execute with step 1 through PROCRUSTES
manopt_sep_start_time = tic();
transf_manopt_gen = som_manopt_genproc(T_globalframe_nois, Tijs_vec_nois, edges, params, matStack(transf_initguess));
% procrustes_end_time = tic();
exectime_manopt_sep = toc(manopt_sep_start_time);

%4) compare output results

testdata.gi = transf_manopt_sep;
% [rotErr,translErr,scale_ratio,translErrNorm] = testNetworkComputeErrors(testdata)
[rotation_error_manopt_sep,translation_error_manopt_sep] = testNetworkComputeErrors(testdata);

testdata.gi = transf_manopt_gen;
[rotation_error_manopt_gen,translation_error_manopt_gen] = testNetworkComputeErrors(testdata);


end %function

