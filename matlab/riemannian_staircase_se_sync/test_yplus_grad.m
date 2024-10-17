clc;
clear;
close all;

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
hessian_mode = 'auto'; 
initguess_is_available = boolean(1);
som_params = struct('N', N, 'd', d, 'd_aff', d_aff, ...
    'global_camera_id', global_camera_id, ...
    'num_tests_per_sigma', num_tests_per_sigma, 'transf_end_thresh', transf_end_thresh, ...
    'max_icp_iterations', max_icp_iterations, 'num_edges_full', num_edges_full, ...
    'num_edges', num_edges, 'procrustes_mode', procrustes_mode, ...
    'riem_grad_mode', riem_grad_mode, ...
    'hessian_mode', hessian_mode, ...
    'initguess_is_available', initguess_is_available);

sigma = 0.0;

edges = (testdata.E);

%set data (no noise)
Tijs_vec = G2T(testdata.gijtruth);
T_globalframe = G2T(testdata.gitruth);

%add noise to data
Tijs_vec_nois = Tijs_vec + sigma.*randn(size(Tijs_vec));
T_globalframe_nois = T_globalframe + sigma.*randn(size(T_globalframe));

R_truth=G2R(testdata.gitruth);
vR_noise=rot_randTangentNormVector(R_truth);
%R_initguess = G2R(rot_randn(testdata.gitruth, sigma_init, N));
R_initguess=rot_exp(R_truth,sigma*pi/5*vR_noise);
transl_initguess = T_globalframe + sigma.*randn(size(T_globalframe));
transf_initguess = RT2G(R_initguess, transl_initguess);
% manopt_start_time = tic();
transf_manopt = som_manopt(T_globalframe_nois, Tijs_vec_nois, edges, som_params, matStack(transf_initguess));
% exectime_manopt = toc(manopt_start_time);


measurements = struct;
measurements.edges = testdata.E;
measurements.R = squeeze(mat2cell(G2R(testdata.gij), som_params.d, som_params.d, ones(som_params.num_edges, 1)))';
measurements.t = mat2cell(G2T(testdata.gij), som_params.d, ones(som_params.num_edges, 1));
measurements.kappa = num2cell(ones(1, som_params.num_edges));
measurements.tau = num2cell(ones(1, som_params.num_edges));
%SoM/testdata fields
measurements.T_globalframe_stiefel = T_globalframe_nois;
measurements.Tijs_vec = Tijs_vec_nois;

[rot_sesync_riemstair, rot_sesync_riemstair_stief, ...
    transl_sesync_riemstair, transl_sesync_riemstair_stiefel, ...
    final_cost_sesync_riemstair, last_num_rows_stiefel] = ...
        som_riemstair_se_sync( ...
            measurements, R_initguess, transl_initguess, som_params);

% figure(1);
% [xt,xDott]=real_randGeodFun(randn(d,d,N));
% fx= @(t) mycost(xt(t),L_T,P_T, cost_const_term_tij);
% dft= @(t) trace( matStack(myeuclgradient(xt(t),L_T,P_T))' * matStack(xDott(t)));
% funCheckDer(fx,dft)

num_rows_stiefel = d;
measurements.d_stiefel = num_rows_stiefel;

problem_data = struct;
problem_data.d = d;

[L_stiefel, P_stiefel] = make_LT_PT_noloops_stiefel( ...
        measurements.T_globalframe_stiefel, measurements.Tijs_vec, measurements.edges, measurements.d_stiefel, som_params);
cost_const_term_tij = compute_fixed_cost_term(measurements.Tijs_vec, problem_data.d);

problem_data = construct_problem_data(measurements);

SE_Sync_opts.r0 = problem_data.d; %Can it actually be d and not d+1?
R_manopt_stiefel = stiefelfactory(SE_Sync_opts.r0,problem_data.d,problem_data.n);
fprintf("R_manopt_stiefel size %g\n", R_manopt_stiefel.dim());
manopt_data.M = R_manopt_stiefel; %M = manifold
myeuclgradient(rot_sesync_riemstair_stief, problem_data, measurements, som_params)
% manopt_data.grad = manopt_data.M.proj(rot_sesync_riemstair_stief, myeuclgradient(rot_sesync_riemstair_stief, problem_data, measurements, som_params))
manopt_data.M.proj(rot_sesync_riemstair_stief, myeuclgradient(rot_sesync_riemstair_stief, problem_data, measurements, som_params))

Yplus = cat_zero_rows_3d_array(rot_sesync_riemstair_stief);
SE_Sync_opts.r0 = problem_data.d; %Can it actually be d and not d+1?
R_manopt_stiefel = stiefelfactory(SE_Sync_opts.r0+1,problem_data.d,problem_data.n);
fprintf("R_manopt_stiefel size %g\n", R_manopt_stiefel.dim());
manopt_data.M = R_manopt_stiefel; %M = manifold
% myeuclgradient(Yplus, problem_data, measurements, som_params)
% manopt_data.grad = manopt_data.M.proj(rot_sesync_riemstair_stief, myeuclgradient(rot_sesync_riemstair_stief, problem_data, measurements, som_params))
Yplus_grad = manopt_data.M.proj(Yplus, myeuclgradient(Yplus, problem_data, measurements, som_params))

eps = 1e-5;
max(Yplus_grad, [], 'all') < eps


%%
function f = mycost(x, problem_data, measurements, som_params)
    measurements.T_globalframe_stiefel = cat_zero_row(measurements.T_globalframe_stiefel, measurements.d_stiefel-problem_data.d);
    [L_stiefel, P_stiefel] = make_LT_PT_noloops_stiefel( ...
        measurements.T_globalframe_stiefel, measurements.Tijs_vec, measurements.edges, measurements.d_stiefel, som_params);
    cost_const_term_tij = compute_fixed_cost_term(measurements.Tijs_vec, problem_data.d);
    f = trace(matStack(x)' * L_stiefel * matStack(x) + matStack(x)' * P_stiefel) + cost_const_term_tij;
end
%%
function g = myeuclgradient(x, problem_data, measurements, som_params)
    measurements.T_globalframe_stiefel = cat_zero_row(measurements.T_globalframe_stiefel, size(x,1)-problem_data.d);
    [L_stiefel, P_stiefel] = make_LT_PT_noloops_stiefel( ...
        measurements.T_globalframe_stiefel, measurements.Tijs_vec, measurements.edges, size(x,1), som_params);
    g = matUnstack(L_stiefel*matStack(x) + (L_stiefel')*matStack(x) + P_stiefel, size(x, 1));
end
%%
% function h = myriemgradient(x, L_stiefel, P_stiefel)
% g = myeuclgradient(x,L_stiefel,P_stiefel);
% h = multiprod(multitransp(x), g) - multiprod(multitransp(g), x);
% h = 0.5 .* h;
% end
%%
function eucl_hess = myeuclhess(x, u, problem_data, measurements, som_params)
measurements.T_globalframe_stiefel = cat_zero_row(measurements.T_globalframe_stiefel, size(x,1)-problem_data.d);
[L_stiefel, P_stiefel] = make_LT_PT_noloops_stiefel( ...
        measurements.T_globalframe_stiefel, measurements.Tijs_vec, measurements.edges, size(x,1), som_params);
P_3d = matUnstack(P_stiefel, size(x, 1));
eucl_hess = multiprod3(u, multitransp(x), P_3d) ...
    + multiprod3(x, multitransp(u), P_3d) ...
    - multiprod3(u, multitransp(P_3d), x) ...
    - multiprod3(x, multitransp(P_3d), u);
eucl_hess = 0.5 .* eucl_hess;
end
%%
% function hess = myhess(x,u,L,P)
% g = myeuclhess(x,u,L,P);
% hess = multiprod(multitransp(x), g) - multiprod(multitransp(g), x);
% hess = 0.5.*hess;
% end
%%
