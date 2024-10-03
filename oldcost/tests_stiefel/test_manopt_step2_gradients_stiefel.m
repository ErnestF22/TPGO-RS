%in order to check whether rotation setup is ok -> check if A*T-b is = 0 

% \begin{equation} \label{eq:coord_desc_second_step}
%     \min_T \norm{A(R) \vec{T} - b(R)}^2,
% \end{equation} 
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




% 1a) PW TRANSLATION DATA INPUT: R, T are the gt, Tijs_mat_nois are the input data
%Rijs, Tijs are the used just for making these dummy tests without having
%the real data

edges = (testdata.E); %graph is not full by default

Tijs_vec = G2T(testdata.gijtruth);
T_globalframe = G2T(testdata.gitruth);

Tijs_mat = tijs_vec_2_tijs_mat(Tijs_vec, edges, N);
%T DATA INPUT (w.r.t. chosen camera id frame)

R_initguess = G2R(testdata.gitruth);
R_initguess_stiefel = cat_zero_rows_3d_array(R_initguess);
transl_initguess = T_globalframe;
transf_gt = testdata.gitruth;

T_globalframe_nois = T_globalframe;
T_globalframe_nois_stiefel = cat_zero_row(T_globalframe_nois);
Tijs_vec_nois = Tijs_vec;
Tijs_mat_nois = Tijs_mat;

% 2a) make L(T), P(T) matrices
% [L_stiefel, P_stiefel] = make_LT_PT_noloops_stiefel( ...
%     T_globalframe_nois_stiefel, Tijs_vec_nois, edges, d+1, som_params);


disp("NOTE: in som_simple_coord_desc() at least one iteration is done by default");
%COORD DESC - step 1
% cost_const_term_tij = compute_fixed_cost_term(Tijs_vec_nois, d);
num_rows_stiefel = d+1;

% Avoid errors when computing R
R = G2R(testdata.gitruth);
R_stiefel = cat_zero_rows_3d_array(R);

% 2b) make A(R), b(R) matrices
[A_stiefel,b] = make_A_b_noloops_stiefel(R_stiefel, T_globalframe_nois, Tijs_vec_nois, edges, num_rows_stiefel, som_params);


%Create a parametric line x(t) in a random direction

% figure(1);
% [xt,xDott]=real_randGeodFun(randn(d,d,N));
% fx= @(t) mycost(xt(t),L_T,P_T, cost_const_term_tij);
% dft= @(t) trace( matStack(myeuclgradient(xt(t),L_T,P_T))' * matStack(xDott(t)));
% funCheckDer(fx,dft)

R_manopt_stiefel = euclideanfactory(num_rows_stiefel,N);
step1.M = R_manopt_stiefel; %M = manifold
% step1.grad = @(x) step1.M.proj(x, myeuclgradient(x, L_stiefel, P_stiefel));

step1.cost = @(x) mycost(x, A_stiefel, b);

figure(1);
title("Transl. gradient check()");

%step1.egrad = @(x) (L_T+L_T')*x + P_T;
step1.egrad = @(x) myeuclgradient(x, A_stiefel, b);
%     step1.grad = @(x) R_manopt_stiefel.egrad2rgrad(matStack(x), step1.egrad);
step1.grad = @(x) step1.M.proj(x, myeuclgradient(x, A_stiefel, b));

checkgradient(step1);

figure(2);
title("Transl. Hessian check()");

step1.hess = @(x, u) step1.M.proj(x, myeuclhess(x, u, A_stiefel, b));

checkhessian(step1);



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

