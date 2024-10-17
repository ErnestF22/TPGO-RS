function test_hessian_manopt_genproc

clc;
clear;
close all;

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
initguess_is_available = boolean(1);
som_params = struct('N', N, 'd', d, 'd_aff', d_aff, ...
    'global_camera_id', global_camera_id, ...
    'num_tests_per_sigma', num_tests_per_sigma, 'transf_end_thresh', transf_end_thresh, ...
    'max_icp_iterations', max_icp_iterations, 'num_edges_full', num_edges_full, ...
    'num_edges', num_edges, 'procrustes_mode', procrustes_mode, ...
    'riem_grad_mode', riem_grad_mode, ...
    'initguess_is_available', initguess_is_available);



% Construct a manifold structure representing the product of groups of
% rotations with the Euclidean space for A. We optimize simultaneously
% for the reference cloud and for the rotations that affect each of the
% measured clouds. Notice that there is a group invariance because
% there is no way of telling which orientation the reference cloud
% should be in.
tuple.R = rotationsfactory(d, N);
tuple.A = euclideanfactory(d, N);
M = productmanifold(tuple);


edges = (testdata.E);
Tijs_vec = G2T(testdata.gijtruth);
T_globalframe = G2T(testdata.gitruth);
Tijs_mat = tijs_vec_2_tijs_mat(Tijs_vec, edges, N);

% Define the cost function here. Points on the manifold M are
% structures with fields X.A and X.R, containing matrices of sizes
% respectively nxm and nxnxN. The store structure (the caching system)
% is used to keep the residue matrix E in memory, as it is also used in
% the computation of the gradient and of the Hessian. This way, we
% prevent redundant computations.
function [f, store] = cost(X, store)
    if ~isfield(store, 'E')
        R = X.R;
        A = X.A;
        store.E = multiprod(R, A);
    end
    %L,P,const_cost_term_tij
    [L, P] = make_LT_PT_noloops(T_globalframe, Tijs_vec, edges, som_params);
    cost_const_term_tij = compute_fixed_cost_term(Tijs_vec, d); 
    f = trace((matStack(R))'*L*(matStack(R)) + (matStack(R))'*P) + cost_const_term_tij;
    %A,b
    %Note: variable name A is already used for the translations!
    %[Amat,b] = make_A_b_noloops(R, T_globalframe, Tijs_vec, Tijs_mat, edges, params);
    %f.A = (Amat*A(:) + b)'*(Amat*A(:) + b);
end

% Riemannian gradient of the cost function.
function [g, store] = grad(X, store)
    R = X.R;
    A = X.A;
    if ~isfield(store, 'E')
        [~, store] = cost(X, store);
    end
    E = store.E;
    %make L,P
    [L, P] = make_LT_PT_noloops(T_globalframe, Tijs_vec, edges, som_params);
    grad_R_1 = (L*matStack(R) + (L')*matStack(R)) + P;
    grad_R_2 = matUnstack(grad_R_1);
    %     h = zeros(size(g));
    g_part = multiprod(multitransp(R),grad_R_2);
    g_asym = (g_part - multitransp(g_part));
    g.R = multiprod(R, g_asym);
    %A,b
    %Note: variable name A is already used for the translations!
    [Amat,b] = make_A_b_noloops(R, T_globalframe, Tijs_vec, Tijs_mat, edges, som_params);
    g.A = reshape(2*(Amat' * Amat) * A(:) + 2*Amat'*b, d, N);
    egrad.R = matUnstack(L*matStack(R) + (L')*matStack(R) + P);
    egrad.A = reshape(2 * (Amat' * Amat) * A(:) + 2*Amat'*b,d,[]);
    store.egrad = egrad;
end

% It is not necessary to define the Hessian of the cost. We do it
% mostly to illustrate how to do it and to study the spectrum of the
% Hessian at the solution (see further down).
function [h, store] = hess(X, d, store)
    R = X.R;
    A = X.A;    

    d_R = tuple.R.tangent2ambient(R, d.R);
    d_A = d.A;
    
    [~, P] = make_LT_PT_noloops(T_globalframe, Tijs_vec, edges, som_params);
    P_3d = matUnstack(P);
    hess = multiprod3(d_R, multitransp(R),P_3d) ...
        - multiprod3(d_R, multitransp(P_3d),R) ...
        - 2.*multiprod3(R, multitransp(P_3d), d_R) ...
        + 2.*multiprod3(R, multitransp(d_R), P_3d) ...
        - multiprod3(multiprod3(R, multitransp(P_3d), R), multitransp(d_R), R) ...
        + multiprod3(P_3d, multitransp(d_R), R);
    const = 0.5;
    hess = const.*hess;
    h.R = hess;

    [Amat,~] = make_A_b_noloops(R, T_globalframe, Tijs_vec, Tijs_mat, edges, som_params);

    h.A = 2.*(Amat'*Amat)*d_A(:);

    store.hess = h;
end


% Setup the problem structure with manifold M and cost+grad functions.
problem.M = M;
problem.cost = @cost;
problem.grad = @grad;
problem.hess = @hess;


% An alternative way to compute the gradient and the hessian is to use 
% automatic differentiation provided in the deep learning toolbox (slower)
% problem.cost = @cost_AD;
%    function f = cost_AD(X)
%        R = X.R;
%        A = X.A;
%        E = multiprod(R, A) - A_measure;
%        f = (E(:)'*E(:))/(2*N);
%    end
% call manoptAD to prepare AD for the problem structure
% problem = manoptAD(problem);


disp("NOTE: At least one ICP iteration is done by default");

% For debugging, it's always nice to check the gradient a few times.
figure(1);
checkgradient(problem);
fprintf("\n");
% pause;
figure(2);
checkhessian(problem);
% pause;

sigma = 0.0;
R_truth = G2R(testdata.gitruth);
vR_noise=rot_randTangentNormVector(R_truth);
R_initguess=rot_exp(R_truth,sigma*pi/5*vR_noise);
transl_initguess = T_globalframe + sigma.*randn(size(T_globalframe));
transf_initguess.R = R_initguess;
transf_initguess.A = transl_initguess;

% Call a solver on our problem. This can probably be much improved if a
% clever initial guess is used instead of a random one.
options.verbosity = 0;
if initguess_is_available
    X = trustregions(problem, transf_initguess, options);
else
    X = trustregions(problem, [], options); %[] is in place of initguess
end
A = X.A;
R = X.R;

transf_prev_mat = testdata.gitruth;
transf_out = RT2G(R, A);


fprintf("\n\n");

end %file function
