function test_step2_gradients_geodesic_stiefel

testdata = testNetwork_som(3); %4 would be the default


% 0a) SOM PARAMS
% som = ShapeOfMotion('noise_test_params.csv'); %params reading is done directly in constructor
%copy the list below from the properties list
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
riem_grad_mode = 'auto'; %'auto' or 'manual'
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

% 0b) SOM PARAMS
%NOTE: sigmas, mus can be seen as couples for each test
sigmas = readmatrix("sigmas.txt"); %sigma = stdev, sigma.^2 = variance
mus = readmatrix("../../data/mus.txt"); %OBS. generally, mus can be d-dimensional; here, we just assume them as scalar (i.e. a d-dimensional vector with all coordinates equal)

%TOUSE: when multple sigmas and mus, and multiple tests per each pair
% for ii = 1:size(sigmas,1)
%     noise_params = struct('sigma', sigmas(ii), 'mu', mus(ii)); 
%     do_som();
% end    
noise_params = struct('sigma', sigmas(1), 'mu', mus(1)); 
sigma = noise_params.sigma;
mu = noise_params.mu;


% 1a) PW TRANSLATION DATA INPUT: R, T are the gt, Tijs_mat_nois are the input data
%Rijs, Tijs are the used just for making these dummy tests without having
%the real data

edges = (testdata.E); %graph is not full by default

Tijs_vec = G2T(testdata.gijtruth);
T_globalframe = G2T(testdata.gitruth);
% Tijs_mat = tijs_vec_2_tijs_mat(Tijs_vec, edges, N);
%T DATA INPUT (w.r.t. chosen camera id frame)

R_initguess = G2R(testdata.gitruth);
% transl_initguess = T_globalframe;


T_globalframe_nois = T_globalframe;
Tijs_vec_nois = Tijs_vec;
% Tijs_mat_nois = Tijs_mat;


% [L_T, P_T] = make_LT_PT_noloops(T_globalframe_nois, Tijs_vec_nois, edges, som_params);


%COORD DESC - step 1
cost_const_term_tij = compute_fixed_cost_term(Tijs_vec_nois, d);
num_rows_stiefel = d;
[R, R_cost, R_info, R_options] = som_stepone_manopt_stiefel(T_globalframe, Tijs_vec, edges, cost_const_term_tij, num_rows_stiefel, som_params, R_initguess);
R_stiefel = cat_zero_rows_3d_array(R);

% Avoid errors when computing R
% R = R_globalframe;
[A_stiefel,b] = make_A_b_noloops_stiefel(R, T_globalframe_nois, Tijs_vec_nois, edges, num_rows_stiefel, som_params);


%COORD DESC - step 2
% [T, T_cost, T_info, T_options] = som_steptwo_manopt(R, T_globalframe_nois, Tijs_vec_nois, edges, som_params, T_initguess);
% [T,A,b] = som_steptwo_procrustes(R, T_globalframe_nois, Tijs_tmp_nois, Tijs_nois, correspondences, som_params);
% disp(T);

% Avoid errors when computing T
T = T_globalframe(:);

%With exact rotation and no noise, the entire vector below should be equal to 0
% A*T_globalframe_nois(:) + b; %== zeros(size(b));

g_T = 2*(A_stiefel' * A_stiefel) * T + 2.*A_stiefel'*b %g = Euclidean gradient


figure(1);
title("Transl. gradient geodesic test");
[xt,xDott]=real_randGeodFun(randn(num_rows_stiefel*N,1));
f_transl= @(t) mycost(xt(t),A_stiefel,b);
dft= @(t) myeuclgradient(xt(t),A_stiefel,b)' * xDott(t);
funCheckDer(f_transl,dft)


figure(2);
title("Transl. Hessian geodesic test");
ddft = @(t) myeuclhess(xt(t), xDott(t), A_stiefel, b);
hess_transl = @(t) xDott(t)' * ddft(t);
funCheckDer(dft, hess_transl)

end %file function

%%
function f = mycost(x, A_stiefel, b)
    e=A_stiefel*x(:) + b;
    f = e'*e;
end
%%
function g = myeuclgradient(x, A_stiefel, b) 
    g = 2*(A_stiefel' * A_stiefel) * x(:) + 2*A_stiefel'*b;
end
%%
function eucl_hess = myeuclhess(x,u,A_stiefel,b)
    eucl_hess = 2.*(A_stiefel'*A_stiefel)*u(:);
end


