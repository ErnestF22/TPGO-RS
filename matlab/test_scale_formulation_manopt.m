function test_scale_formulation_manopt

% importmanopt;


%% testNetwork

N = 5;
mindeg = 3;
testdata = testNetwork_params(3, N, 'banded', mindeg);
testdata.mindeg = 3;

X_gt.lambda = testdata.lambdaijtruth;
X_gt.R = G2R(testdata.gitruth);
X_gt.T = G2T(testdata.gitruth);

disp("X_gt")
disp(X_gt)

tijs = G2T(testdata.gij);
edges = testdata.E;


%%

nrs = size(X_gt.T, 1);
d = size(tijs, 1); % d = 3
N = size(X_gt.R, 3);

sz=[nrs,d,N];

num_edges = size(edges, 1);


rho = 1.0; %TODO: set this properly

% variables random generation/init
tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
tuple.lambda = euclideanfactory(num_edges, 1);
M = productmanifold(tuple);

problem_data = struct('sz', sz, 'edges', edges, 'tijs', tijs);

problem_data.rho = rho; 

%% output (problem_curve_data) definition


problem.M = M;
problem.cost = @(x) ssom_cost(x, problem_data);
% problem.egrad = @(x) ssom_egrad(x, problem_data);
% problem.ehess = @(x, u) ssom_ehess_genproc(x, u, problem_data);
problem.grad = @(x) ssom_rgrad(x, problem_data);
problem.hess = @(x, u) ssom_rhess_genproc(x, u, problem_data);

fprintf("\n");
disp("Gradient check easy");
X_gradcheck.R = eye3d(nrs, d, N);
X_gradcheck.T = ones(nrs, N);
X_gradcheck.lambda = zeros(num_edges, 1);
figure(1)
checkgradient(problem, X_gradcheck);

fprintf("\n");
disp("Gradient check rand");
figure(2)
checkgradient(problem);

close all;

fprintf("\n");
disp("Hessian check easy");
% X_gradcheck.R = eye3d(nrs, d, N);
% X_gradcheck.T = zeros(nrs, N);
% X_gradcheck.lambda = zeros(num_edges, 1);
% Xdot_gradcheck.R = eye3d(nrs, d, N);
% Xdot_gradcheck.T = zeros(nrs, N);
% Xdot_gradcheck.lambda = zeros(num_edges, 1);
figure(3)
X_hesscheck = M.rand();
Xdot_hesscheck = M.rand();
X_hesscheck.R = eye3d(nrs, d, N);
Xdot_hesscheck.R = stiefel_randTangentNormVector(X_hesscheck.R);
% only case this works is with the following T, lambda and R = zeros()
X_hesscheck.T = ones(nrs, N);
X_hesscheck.lambda = 2 * ones(num_edges, 1);
checkhessian(problem, X_hesscheck, Xdot_hesscheck);


fprintf("\n");
disp("Hessian check rand");
checkhessian(problem);

% problem_data.R_gt = eye3d(nrs, d, N);
% problem_data.T_gt = zeros(nrs, N);
% problem_data.lambda_gt = zeros(num_edges, 1);



% trustregions(problem);


end %file function

%%
function [h] = ssom_rhess_genproc(X, Xdot, problem_data)
R = X.R;
T = X.T;
lambdas = X.lambda;
Rdot = Xdot.R;
Tdot = Xdot.T;
lambdasdot = Xdot.lambda;

% Careful: tangent vectors on the rotation group are represented as
% skew symmetric matrices. To obtain the corresponding vectors in
% the ambient space, we need a little transformation. This
% transformation is typically not needed when we compute the
% formulas for the gradient and the Hessian directly in Riemannian
% form instead of resorting the egrad2rgrad and ehess2rhess. These
% latter tools are convenient for prototyping but are not always
% the most efficient form to execute the computations.

hrt = ssom_ehess_r_t(R, T, Tdot, lambdas, problem_data);
htr = ssom_ehess_t_r(R, Rdot, T, lambdas, problem_data);

% h_lambda_lambda = zeros(size(lambda));
h_lambda_lambda = ssom_ehess_lambda_lambda(lambdas, lambdasdot, R, problem_data);

% h_r_lambda = zeros(size(hrt));
h_r_lambda = ssom_ehess_r_lambda(R, T, lambdas, lambdasdot, problem_data);

% h_t_lambda = zeros(size(htr));
h_t_lambda = ssom_ehess_t_lambda(R, T, lambdas, lambdasdot, problem_data);

% h_lambda_r = zeros(size(h_lambda_lambda));
h_lambda_r = ssom_ehess_lambda_r(X, Xdot, problem_data);

% h_lambda_t = zeros(size(h_lambda_lambda));
h_lambda_t = ssom_ehess_lambda_t(R, T, Tdot, lambdas, problem_data);

ehR = ssom_ehess_R_R(R, Rdot, problem_data) + hrt + h_r_lambda;
egR = egrad_R(R, T, lambdas, problem_data);
h.R = manopt_stiefel_ehess2rhess(R, egR, ehR, Rdot);
h.T = ssom_ehess_T_T(R, T, Tdot, lambdas, problem_data) + htr + h_t_lambda;
h.lambda = h_lambda_lambda + h_lambda_r + h_lambda_t;
end %ehess genproc

function h = ssom_ehess_lambda_lambda(x, xdot, R, problem_data)
% h_lambda_lambda = zeros(1,1)
% lambdas_dot = Xdot.lambda;
edges = problem_data.edges;
tijs = problem_data.tijs;
% rho = problem_data.rho;

h = zeros(length(x), 1);

num_edges = size(edges, 1);
for ee = 1:num_edges
    ii = edges(ee, 1);
    % jj = edges(ee, 2);
    % lambda_e = lambdas(ee);
    tij_e = tijs(:, ee);
    % T_i = problem_data.T(:, ii);
    % T_j = problem_data.T(:, jj);
    R_i = R(:, :, ii);
    % a = T_i - T_j;
    b = R_i * tij_e;
    h(ee) = 2*(b' * b) * xdot(ee);
end

end

function h = ssom_ehess_T_T(R, ~, Tdot, lambdas, problem_data)
tijs_scaled = make_tijs_scaled(lambdas, problem_data.tijs);
[LR] = make_LR_PR_BR_noloops(R, tijs_scaled, problem_data.edges);
% gT = T * (problem.LR+problem.LR')+(problem.PR)';
h = Tdot*(LR' + LR);
end

function h = ssom_ehess_R_R(R, ~, ~)
h = zeros(size(R));
end

%% 2
function h = ssom_ehess_r_lambda(~, T, ~, lambdas_dot, problem_data)

nrs = size(T, 1);
d = size(problem_data.tijs, 1);
N = size(T, 2);

Ph = zeros(nrs, d*N);

tijs_dot_scaled = make_tijs_scaled(lambdas_dot, problem_data.tijs);

idx_col_p = reshape(1:d*N, [], N)';

num_edges = size(problem_data.edges,1);
for e = 1:num_edges
    ii = problem_data.edges(e,1);
    jj = problem_data.edges(e,2);
    T_j = T(:, jj);
    T_i = T(:, ii);
    tij_dot = tijs_dot_scaled(:,e);
    P_e = 2 * (T_i * tij_dot' - T_j * tij_dot');
    Ph(:, idx_col_p(ii, :)) = ...
        Ph(:, idx_col_p(ii, :)) + P_e;
end

h=matUnstackH(Ph,d);

end


%% 3
function h = ssom_ehess_t_lambda(R, T, ~, lambdas_dot, problem_data)
% h_t_lambda = zeros(size(htr));

edges = problem_data.edges;
% rho = problem_data.rho;

h = zeros(size(T));
N = size(R,3);
num_edges = size(edges, 1);
for e = 1:num_edges
    ii = edges(e,1);
    jj = edges(e,2);
    Ri = R(:,:,ii);
    %         Rj = X.R(:,:,jj);
    %
    BIJ = zeros(N,1);
    BIJ(ii) = 1;
    BIJ(jj) = -1;
    %
    tij = problem_data.tijs(:, e);
    w_ij = BIJ * lambdas_dot(e) * tij' * Ri';
    h = h + 2 * w_ij';
end

end



%% 4
function eh = ssom_ehess_lambda_r(X, Xdot, problem_data)

R = X.R;
Rdot = Xdot.R;
T = X.T;
lambdas = X.lambda;
% h_lambda_t = zeros(size(h_lambda_lambda));

% x = X.lambda;
% lambdas_dot = Xdot.lambda;
edges = problem_data.edges;
tijs_vec = problem_data.tijs;
% rho = problem_data.rho;

% nrs = problem_data.sz(1);
% % d = problem_data.sz(2);
% N = problem_data.sz(3);

eh = zeros(size(lambdas));

num_edges = size(edges, 1);
for ee = 1:num_edges
    ii = edges(ee, 1);
    jj = edges(ee, 2);
    lambda_e = lambdas(ee);
    tij = tijs_vec(:, ee);
    T_i = T(:, ii);
    T_j = T(:, jj);
    R_i_dot = Rdot(:, :, ii);
    R_i = R(:, :, ii);
    a = T_i - T_j;
    bdot = R_i_dot * tij;
    b = R_i * tij;
    e_th_elem_half = lambda_e * (bdot' * b) + lambda_e * (b' * bdot)  + ...
        a' * bdot;

    eh(ee) = 2 * e_th_elem_half;
end

end


%% 5
function h = ssom_ehess_lambda_t(R, ~,  Tdot, lambdas, problem_data)

% h_lambda_r = zeros(size(h_lambda_lambda));

% x = X.lambdas;
% lambdas_dot = Xdot.lambdas;
edges = problem_data.edges;
tijs_vec = problem_data.tijs;
% rho = problem_data.rho;

h = zeros(size(lambdas));

num_edges = size(edges, 1);
for ee = 1:num_edges
    ii = edges(ee, 1);
    jj = edges(ee, 2);
    % lambda_e = x(ee);
    tij = tijs_vec(:, ee);
    % T_i = T(:,ii);
    % T_j = T(:,jj);
    T_i_dot = Tdot(:, ii);
    T_j_dot = Tdot(:, jj);
    % a = T_i - T_j;
    R_i = R(:, :, ii);
    b = R_i * tij;
    adot = T_i_dot - T_j_dot;
    e_th_elem = 2 * adot' * b;
    h(ee) = e_th_elem;
end

end

function h = ssom_ehess_r_t(~, T, Tdot, lambdas, problem_data)
nrs = size(T, 1);
d = size(problem_data.tijs, 1);
N = size(T, 2);

Ph = zeros(nrs, d*N);

tijs_scaled = make_tijs_scaled(lambdas, problem_data.tijs);

idx_col_p = reshape(1:d*N, [], N)';

num_edges = size(problem_data.edges,1);
for e = 1:num_edges
    ii = problem_data.edges(e,1);
    jj = problem_data.edges(e,2);
    Tj_dot = Tdot(:, jj);
    Ti_dot = Tdot(:, ii);
    tij = tijs_scaled(:,e);
    P_e = 2 * (Ti_dot * tij' - Tj_dot * tij');
    Ph(:, idx_col_p(ii, :)) = ...
        Ph(:, idx_col_p(ii, :)) + P_e;
end

h=matUnstackH(Ph,d);

end

function h = ssom_ehess_t_r(~, dR, ~, lambdas, problem_data)
tijs_scaled = make_tijs_scaled(lambdas, problem_data.tijs);
[~, PR_dot] = make_LR_PR_BR_noloops(dR, tijs_scaled, problem_data.edges);
h=PR_dot';
end

%% RHESS CONVERSION

function rhess = manopt_stiefel_ehess2rhess(X, egrad, ehess, H)
XtG = multiprod(multitransp(X), egrad);
symXtG = multisym(XtG);
HsymXtG = multiprod(H, symXtG);
rhess = stiefel_tangentProj(X, ehess - HsymXtG);
end

function g = egrad_R(~, T, lambdas, problem_data)
nrs = size(T, 1);
d = size(problem_data.tijs, 1);
N = size(T, 2);

P = zeros(nrs, d*N);

tijs_scaled = make_tijs_scaled(lambdas, problem_data.tijs);

idx_col_p = reshape(1:d*N, [], N)';

num_edges = size(problem_data.edges,1);
for e = 1:num_edges
    ii = problem_data.edges(e,1);
    jj = problem_data.edges(e,2); 
    T_j = T(:, jj);
    T_i = T(:, ii);
    tij = tijs_scaled(:,e);
    P_e = 2 * (T_i * tij' - T_j * tij');
    P(:, idx_col_p(ii, :)) = ...
        P(:, idx_col_p(ii, :)) + P_e;
end

g=matUnstackH(P,d);
end