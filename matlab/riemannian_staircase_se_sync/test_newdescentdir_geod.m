close all;
clear;
clc;


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

SE_Sync_opts = struct;

SE_Sync_opts.r0 = d;
if isfield(SE_Sync_opts, 'eig_comp_max_iters')
    fprintf(' Maximum number of iterations to perform for minimum eigenvalue computation in test for positive semidefiniteness: %g\n', SE_Sync_opts.eig_comp_max_iters);
else
    SE_Sync_opts.eig_comp_max_iters = 2000;
    fprintf(' Maximum number of iterations to perform for minimum eigenvalue computation in test for positive semidefiniteness: %g [default]\n', SE_Sync_opts.eig_comp_max_iters);
end
if isfield(SE_Sync_opts, 'rmax')
    fprintf(' Final level of Riemannian Staircase: %d\n', SE_Sync_opts.rmax);
else
    SE_Sync_opts.rmax = 7;
    fprintf(' Setting final level of Riemannian Staircase to %d [default]\n', SE_Sync_opts.rmax);
end
if ~isfield(SE_Sync_opts, 'Cholesky')
    fprintf(' Using QR decomposition to compute orthogonal projection [default]\n');
    SE_Sync_opts.Cholesky = false;
else
    if SE_Sync_opts.Cholesky
        fprintf(' Using Cholesky decomposition to compute orthogonal projection\n');
    else
        fprintf(' Using QR decomposition to compute orthogonal projection\n');
    end
end
if isfield(SE_Sync_opts, 'min_eig_num_tol')
    fprintf(' Tolerance for accepting an eigenvalue as numerically nonnegative in optimality verification: %g\n', SE_Sync_opts.min_eig_num_tol);
else
    SE_Sync_opts.min_eig_num_tol = 1e-4;
    fprintf(' Tolerance for accepting an eigenvalue as numerically nonnegative in optimality verification: %g [default]\n', SE_Sync_opts.min_eig_num_tol);
end

Manopt_opts = struct;
if isfield(Manopt_opts, 'tolgradnorm')
    fprintf(' Stopping tolerance for norm of Riemannian gradient: %g\n', Manopt_opts.tolgradnorm);
else
    Manopt_opts.tolgradnorm = 1e-2;
    fprintf(' Setting stopping tolerance for norm of Riemannian gradient to: %g [default]\n', Manopt_opts.tolgradnorm);
end


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

% [rot_sesync_riemstair, rot_sesync_riemstair_stief, ...
%     transl_sesync_riemstair, transl_sesync_riemstair_stiefel, ...
%     final_cost_sesync_riemstair, last_num_rows_stiefel] = ...
%         som_riemstair_se_sync( ...
%             measurements, R_initguess, transl_initguess, som_params);
% Yopt = rot_sesync_riemstair_stief;

problem_data = struct;
problem_data.d = d;
problem_data.n = N;

R_manopt_stiefel = stiefelfactory(SE_Sync_opts.r0,problem_data.d,problem_data.n);
manopt_data.M = R_manopt_stiefel; %M = manifold
% Define the cost function and its gradient.
% manopt_data.cost = @(x) mycost(x, L_stiefel, P_stiefel, cost_const_term_tij);
measurements.d_stiefel = SE_Sync_opts.r0;
manopt_data.cost = @(x) mycost(x, problem_data, measurements, som_params);
manopt_data.egrad = @(x) myeuclgradient(x, problem_data, measurements, som_params);
manopt_data.grad = @(x) manopt_data.M.proj(x, myeuclgradient(x, problem_data, measurements, som_params));


[Yopt, Fval, manopt_info, Manopt_opts] = manoptsolve(manopt_data, R_initguess, Manopt_opts);
SDPLRval = Fval(end);

num_rows_stiefel = d;
measurements.d_stiefel = num_rows_stiefel;

problem_data = struct;
problem_data.d = d;
problem_data = construct_problem_data(measurements);

%%
figure("Name", "Yopt");
% [xt,xDott]=real_randGeodFun(randn(d,d,N));
% fx= @(t) mycost(xt(t),L_T,P_T, cost_const_term_tij);
% dft= @(t) trace( matStack(myeuclgradient(xt(t),L_T,P_T))' * matStack(xDott(t)));
% funCheckDer(fx,dft)

[stief_t,vt] = stiefel_rand_geod_fun_new(matStack(Yopt));
stief_0 = stief_t(0);
% f_rot= @(t) mycost(stief_t(t),problem_data, measurements, som_params);
% dft_rot= @(t) stiefel_metric
% fplot(@(t) sum(stief_t(t),'all'), [0,1])
% fplot(@(t) sum(stief_0 + t*vt(t),'all'), [0,1])


stief_exp_t = @(t) stiefel_exp(stief_t(t), vt(t));
% Y = @(t) stief_0 + stief_exp_t(t);

fcost = @(t) mycost(stief_t(t), problem_data, measurements, som_params);

%discrete plot
% discrete_interval = (-0.1:0.01:0.1)';
% fcost_discrete_interval = zeros(size(discrete_interval));
% for ii = 1:size(discrete_interval)
%     fcost_discrete_interval(ii) = fcost(discrete_interval(ii));
% end
% plot(discrete_interval, fcost_discrete_interval);
% grid

fplot(@(t) fcost(t), [-0.1,0.1]);
gcost= @(t) matStack(myeuclgradient(matUnstack(stief_t(t), 3), problem_data, measurements, som_params));
dfcost = @(t) stiefel_metric(stief_t(t),gcost(t), vt(t));
figure(1)
funCheckDer(stief_t,vt)
figure(2)
funCheckDer(fcost, dfcost)
keyboard

%% Yplus
Yplus = cat_zero_rows_3d_array(Yopt);
%Computing escape direction
rNext = num_rows_stiefel + 1;
Lambda = compute_Lambda(matStackH(Yopt), problem_data, SE_Sync_opts.Cholesky);

% Compute minimum eigenvalue/eigenvector pair for Q - Lambda
% tic();
[lambda_min, v] = Q_minus_Lambda_min_eig(Lambda, problem_data, Yopt, SE_Sync_opts.min_eig_num_tol, SE_Sync_opts.eig_comp_max_iters, SE_Sync_opts.Cholesky);
% min_eig_comp_time = toc();
disp('Computing escape direction...');
v_stiefel = zeros((rNext)*problem_data.n, 1);
idx_stief = reshape(1:(rNext)*problem_data.n, rNext, problem_data.n);
idx = reshape(1:problem_data.d*problem_data.n, problem_data.d, problem_data.n);
for ii = 1:problem_data.n
    v_stiefel(idx_stief(:,ii)) = vertcat(v(idx(:,ii)), zeros(rNext - problem_data.d, 1));
end
Ydot = horzcat(zeros((rNext) * problem_data.n, problem_data.d-1), v_stiefel);


% Augment the dimensionality of the Stiefel manifolds in
% preparation for the next iteration
%         manopt_data.M = stiefelstackedfactory(problem_data.n, problem_data.d, r+1);
manopt_data.M = stiefelfactory(rNext,problem_data.d,problem_data.n);
measurements.d_stiefel = rNext;
% Update preconditioning function, if it's used
if(exist('precon', 'var'))
    manopt_data.precon = @(x,u) manopt_data.M.proj(x, precon(u));
end

% Perform line search along the escape direction Ydot to escape the
% saddle point and obtain the initial iterate for the next level in
% the Staircase

% Compute a scaling factor alpha such that the scaled step
% alpha*Ydot' should produce a trial point Ytest whose gradient has
% a norm 100 times greater than the gradient tolerance stopping
% criterion currently being used in the RTR optimization routine
alpha = Manopt_opts.tolgradnorm / (norm(v) * abs(lambda_min));

disp('Line searching along escape direction to escape saddle point...');
tic();
%         [stepsize, Y0T] = linesearch_decrease(manopt_data, Yplus', alpha * Ydot', SDPLRval);
manopt_data.cost = @(x) mycost(x, problem_data, measurements, som_params);
manopt_data.egrad = @(x) myeuclgradient(x, problem_data, measurements, som_params);
manopt_data.grad = @(x) manopt_data.M.proj(x, myeuclgradient(x, problem_data, measurements, som_params));
[stepsize, Y02d] = linesearch_decrease(manopt_data, matStack(Yplus), alpha * Ydot, SDPLRval);
line_search_time = toc();
%         Y0 = Y0T';
Y0 = matUnstack(Y02d, rNext);
fprintf('Line search completed (elapsed computation time %g seconds)\n', line_search_time);

figure("Name", "Yplus");
[stief_t_new,vt_new] = ...
    stiefel_rand_geod_fun_new(matStack(Yplus)); %must be initialized with rotation mat

fcost_new = @(t) mycost(stief_t_new(t), problem_data, measurements, som_params);
fplot(@(t) fcost_new(t), [-0.1,0.1]);

figure("Name", "Yplus_new_descent_dir");
stief_exp_t = @(t) stiefel_exp(stief_t_new(t), vt_new(t));
Y = @(t) mycost(matStack(Y0) + stief_exp_t(t), problem_data, measurements, som_params);

% fplot(@(t) Y(t), [-0.1,0.1]);


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
