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

num_rows_stiefel = d+1;

R_initguess = G2R(testdata.gitruth);
R_initguess_stiefel = cat_zero_rows_3d_array(R_initguess);
transl_initguess = T_globalframe;
transf_gt = testdata.gitruth;

T_globalframe_nois = T_globalframe;
T_globalframe_nois_stiefel = cat_zero_row(T_globalframe_nois);
Tijs_vec_nois = Tijs_vec;
Tijs_mat_nois = Tijs_mat;


% 2a) make L(T), P(T) matrices
[L_T_stiefel, P_T_stiefel] = make_LT_PT_noloops_stiefel( ...
    T_globalframe_nois_stiefel, Tijs_vec_nois, edges, num_rows_stiefel, som_params);


disp("NOTE: in som_simple_coord_desc() at least one iteration is done by default");
%COORD DESC - step 1
cost_const_term_tij = compute_fixed_cost_term(Tijs_vec_nois, d);
[R, R_cost, R_info, R_options] = som_stepone_manopt_stiefel(T_globalframe, Tijs_vec, edges, cost_const_term_tij, num_rows_stiefel, som_params, R_initguess_stiefel);
% disp(R);

% Avoid errors when computing R
R = G2R(testdata.gitruth);
R_stiefel = cat_zero_rows_3d_array(R);


g_R = L_T_stiefel*matStack(R_stiefel) + (L_T_stiefel')*matStack(R_stiefel) + P_T_stiefel %g = Euclidean gradient (stacked)
h_R = zeros(size(g_R)); %h = Riemannian gradient
for ii = 1:N
    gpart = g_R(ii*d-(d-1):ii*d, :);
    g_asym = 0.5 .* (gpart - gpart');
    h_R(ii*d-(d-1):ii*d, :) = g_asym';
end
h_R %h = Riemannian gradient

%Create a parametric line x(t) in a random direction

% figure(1);
% [xt,xDott]=real_randGeodFun(randn(d,d,N));
% fx= @(t) mycost(xt(t),L_T,P_T, cost_const_term_tij);
% dft= @(t) trace( matStack(myeuclgradient(xt(t),L_T,P_T))' * matStack(xDott(t)));
% funCheckDer(fx,dft)

R_manopt_stiefel = stiefelfactory(d+1*d,N);
step1.M = R_manopt_stiefel; %M = manifold
% step1.grad = @(x) step1.M.proj(x, myeuclgradient(x, L_stiefel, P_stiefel));

% figure(1);
[stief_t,stief_Dott] = stiefel_randGeodFun(matStack(make_rand_stiefel_mat(d, d, N))); %must be initialized with rotation mat
f_stief = @(t) mycost(matUnstack(stief_t(t)), ...
    L_T_stiefel, P_T_stiefel, cost_const_term_tij);
gEucl_stief = @(t) myeuclgradient(stief_t(t), L_T_stiefel, P_T_stiefel);
g_stief= @(t) step1.M.proj(cat_zero_rows_3d_array(matUnstack(stief_t(t), 3)), gEucl_stief(t));
dft_stief = @(t) sum(stiefel_metric( cat_zero_rows_3d_array(matUnstack(stief_t(t),3)), ...
    gEucl_stief(t), ...
    cat_zero_rows_3d_array(matUnstack(stief_Dott(t), 3)) ));
title("Rot. gradient geod. check");
funCheckDer(f_stief,dft_stief)


% figure(2);
% title("Rot. Hessian geod. check");
% ddft_rot= @(t) myhess(...
%     cat_zero_rows_3d_array(rot_t(t)), ...
%     cat_zero_rows_3d_array(rot_Dott(t)), L_T_stiefel, P_T_stiefel);
% hess_rot= @(t) stiefel_metric(matStack(ddft_rot(t)), ...
%     matStack(cat_zero_rows_3d_array(rot_Dott(t))), ...
%     matStack(ddft_rot(t)));
% funCheckDer(dft_rot, hess_rot)




%%
function f = mycost(x, L, P, fixed_cost_term)
    x = cat_zero_rows_3d_array(x);
    f = trace((matStack(x))'*L*(matStack(x)) + (matStack(x))'*P) + fixed_cost_term;
%     f = trace((matStack(x))'*P) + fixed_cost_term;
end
%%
function g = myeuclgradient(x, L, P) 
    x = matUnstack(x, 3);
    x = cat_zero_rows_3d_array(x);
    g = matUnstack(L*matStack(x) + (L')*matStack(x) + P, size(x, 1));
%     g = P;
end
%%
% function h = myriemgradient(x, L, P) 
%     %g = myeuclgradient(x,L,P) BUT still stacked
%     g = (L*matStack(x) + (L')*matStack(x)) + P;
%     g = matUnstack(g);
% %     h = zeros(size(g));
%     g_part = multiprod(multitransp(x),g);
%     g_asym = (g_part - multitransp(g_part));
%     h = multiprod(x, g_asym);
% end
%%
function eucl_hess = myeuclhess(x, u, L, P)
x = matUnstack(x, 3);
x = cat_zero_rows_3d_array(x);
P_3d = matUnstack(P, size(x, 1));
eucl_hess = multiprod3(u, multitransp(x), P_3d) ... 
        + multiprod3(x, multitransp(u), P_3d) ...
        - multiprod3(u, multitransp(P_3d), x) ...
        - multiprod3(x, multitransp(P_3d), u);
eucl_hess = 0.5 .* eucl_hess;
end
%%
function hess = myhess(x,u,L,P)
g = myeuclhess(x,u,L,P);
hess = multiprod(multitransp(x), g) - multiprod(multitransp(g), x);
hess = 0.5.*hess;
end
%%

% [rot_t,rot_Dott] = rot_randGeodFun(rot_randn(repmat(eye(d),1,1,N))); %must be initialized with rotation mat
% f_rot= @(t) mycost(cat_zero_rows_3d_array(rot_t(t)), ...
%     L_T_stiefel, P_T_stiefel, cost_const_term_tij);
% dft_rot= @(t) stiefel_metric(matStack(cat_zero_rows_3d_array(rot_t(t))), ...
%     matStack(cat_zero_rows_3d_array(rot_Dott(t))), ...
%     stiefel_eucl2RiemGrad(matStack(cat_zero_rows_3d_array(rot_t(t))), myeuclgradient(matStack(cat_zero_rows_3d_array(rot_t(t))), L_T_stiefel, P_T_stiefel)) );
% funCheckDer(f_rot,dft_rot)