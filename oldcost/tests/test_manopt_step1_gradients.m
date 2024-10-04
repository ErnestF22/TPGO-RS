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
riem_grad_mode = 'manual'; %'auto' or 'manual'
hessian_mode = 'manual'; %'auto' or 'manual'
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

Tijs_mat = tijs_vec_2_tijs_mat(Tijs_vec, edges, N);
%T DATA INPUT (w.r.t. chosen camera id frame)

R_initguess = G2R(testdata.gitruth);
transl_initguess = T_globalframe;
transf_initguess.R = R_initguess;
transf_initguess.A = transl_initguess;
transf_gt = testdata.gitruth;

T_globalframe_nois = T_globalframe;
Tijs_vec_nois = Tijs_vec;
Tijs_mat_nois = Tijs_mat;

% 2a) make L(T), P(T) matrices
[L_T, P_T] = make_LT_PT_noloops(T_globalframe_nois, Tijs_vec_nois, edges, som_params);


disp("NOTE: in som_simple_coord_desc() at least one iteration is done by default");
%COORD DESC - step 1
cost_const_term_tij = compute_fixed_cost_term(Tijs_vec_nois, d);
[R, R_cost, R_info, R_options] = som_stepone_manopt(T_globalframe, Tijs_vec, edges, cost_const_term_tij, som_params, R_initguess);
% R = som_stepone_procrustes(T_globalframe_nois, Tijs_tmp_nois, Tijs_nois, correspondences, som_params);
% disp(R);

% Avoid errors when computing R
R = G2R(testdata.gitruth);

[A,b,R] = make_A_b_noloops(R, T_globalframe_nois, Tijs_vec_nois, edges, som_params);

g_R = L_T*matStack(R) + (L_T')*matStack(R) + P_T %g = Euclidean gradient (stacked)
h_R = zeros(size(g_R)); %h = Riemannian gradient
for ii = 1:N
    gpart = g_R(ii*d-(d-1):ii*d, :);
    g_asym = 0.5 .* (gpart - gpart');
    h_R(ii*d-(d-1):ii*d, :) = g_asym';
end
h_R %h = Riemannian gradient

%Create a parametric line x(t) in a random direction

figure(1);
[xt,xDott]=real_randGeodFun(randn(d,d,N));
fx= @(t) mycost(xt(t),L_T,P_T, cost_const_term_tij);
dft= @(t) trace( matStack(myeuclgradient(xt(t),L_T,P_T))' * matStack(xDott(t)));
funCheckDer(fx,dft)

figure(2);
[rot_t,rot_Dott] = rot_randGeodFun(rot_randn(repmat(eye(d),1,1,N))); %must be initialized with rotation mat
f_rot= @(t) mycost(rot_t(t),L_T,P_T, cost_const_term_tij);
dft_rot= @(t) sum(rot_metric(rot_t(t), rot_Dott(t), myriemgradient(rot_t(t), L_T, P_T)));
funCheckDer(f_rot,dft_rot)


% R = rot_randn(repmat(eye(d),1,1,N));
% vrand=rot_randTangentNormVector(R);
% rot_metric(R,vrand,myriemgradient(R, L_T, P_T))
% rot_metric(R,vrand,myeuclgradient(R, L_T, P_T))


%COORD DESC - step 2
% [T, T_cost, T_info, T_options] = som_steptwo_manopt(R, T_globalframe_nois, Tijs_tmp_nois, Tijs_nois, correspondences, som_params);
% % [T,A,b] = som_steptwo_procrustes(R, T_globalframe_nois, Tijs_tmp_nois, Tijs_nois, correspondences, som_params);
% % disp(T);
% 
% 
% transf_curr = make_transf(R, T);
% % disp(transf_curr);

%%
function f = mycost(x, L, P, fixed_cost_term) 
    f = trace((matStack(x))'*L*(matStack(x)) + (matStack(x))'*P) + fixed_cost_term;
end
%%
function g = myeuclgradient(x, L, P) 
    g = matUnstack(L*matStack(x) + (L')*matStack(x) + P);
end
%%
function h = myriemgradient(x, L, P) 
    %g = myeuclgradient(x,L,P) BUT still stacked
    g = (L*matStack(x) + (L')*matStack(x)) + P;
    g = matUnstack(g);
%     h = zeros(size(g));
    g_part = multiprod(multitransp(x),g);
    g_asym = (g_part - multitransp(g_part));
    h = multiprod(x, g_asym);
end
%%

