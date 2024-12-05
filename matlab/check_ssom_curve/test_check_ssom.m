function problem_data=test_check_ssom()

nrs = 3;
d = 3;
N = 5;

sz=[nrs,d,N];

%graph random init
num_edges = 8;
G = graph(true(N), 'omitselfloops'); % Alternative without self-loops
p = randperm(numedges(G), num_edges);
G = graph(G.Edges(p, :));
edges = table2array(G.Edges);

tijs = 10 * rand(d, num_edges);

lambdas = 10 * rand(num_edges, 1);

rho = 1.0; %TODO: make this rand() later

tijs_scaled = make_tijs_scaled(lambdas, tijs); %!!
problem_data = struct('sz', sz, 'edges', edges, 'tijs', tijs, 'Tijs', tijs_scaled);

% variables random generation/init
tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
tuple.lambda = euclideanfactory(num_edges, 1);
M = productmanifold(tuple);
% problem.M = M;
x = M.rand();
% lambdas = x.lambda;
T = x.T;
R = x.R;

%g.R
[problem_data.P, problem_data.frct] = ...
    make_step1_p_fct(T, tijs_scaled, edges);
%g.T
problem_data.Tijs = tijs_scaled;
[problem_data.LR, problem_data.PR, problem_data.BR] = ...
    make_LR_PR_BR_noloops(R, tijs_scaled, edges);
%g.lambda
problem_data.T = T;
problem_data.R = R;
[aL, bL, cL] = makeABClambda(x, problem_data);
problem_data.aL = aL; problem_data.bL = bL; problem_data.cL = cL;

problem_data.rho = rho; %ReLU() part should not be needed for Hessian tests

%% output (problem_curve_data) definition

problem_data.edges = edges;
problem_data.Tijs = tijs;
problem_data.cost_lambda=@(x) ssom_cost_lambda(x,problem_data);
problem_data.cost_R=@(x) rsom_cost_rot_stiefel(x,problem_data);
problem_data.cost_T=@(x) rsom_cost_transl_stiefel(x,problem_data);
%
% problem_data.egrad_lambda=@(x) egrad_lambda(x,problem_data);
problem_data.grad_lambda=@(lambdas,R,T) grad_lambda(lambdas,R,T,problem_data);
problem_data.rgrad_R=@(x) rgrad_R(x,problem_data);
problem_data.egrad_T=@(x) egrad_T(x,problem_data);
%
problem_data.ssom_ehess_lambda_lambda=@(x,u,R) ssom_ehess_lambda_lambda(x,u,R,problem_data);
problem_data.ssom_ehess_r_lambda=@(R, Rdot, T, lambdadot) ssom_ehess_r_lambda(R, Rdot, T, lambdadot, problem_data);
problem_data.ssom_ehess_t_lambda=@(x,u,R,lambdas_dot) ssom_ehess_t_lambda(x,u,R,lambdas_dot,problem_data);
problem_data.ssom_ehess_lambda_r=@(x,u,T,R,Rdot) ssom_ehess_lambda_r(x, u, T, R, Rdot, problem_data);
problem_data.ssom_ehess_lambda_t=@(x,u,T,Tdot,R) ssom_ehess_lambda_t(x,u,T,Tdot,R,problem_data);

end %file function

function g=rgrad_R(x, problem_data)
g = rsom_rgrad_rot_stiefel(x, problem_data);
end

function g=egrad_T(x,problem_data)
g = rsom_egrad_transl_stiefel(x, problem_data);
end

%% 0) grad lambda
% function g=egrad_lambda(x,problem_data)
% %Note: already checked (genproc-wise) through Manopt
% 
% % x in this context is already lambda
% 
% % g = ssom_grad_lambda(x, problem_data);
% 
% edges = problem_data.edges;
% rho = problem_data.rho;
% 
% % x = lambdas in this context
% 
% g = zeros(length(x), 1);
% 
% num_edges = size(edges, 1);
% for ee = 1:num_edges
%     % ii = edges(ee, 1);
%     % jj = edges(ee, 2);
%     lambda_e = x(ee);
%     % tij_e = Tijs_vec(:, ee);
%     % aLi = problem_data.aL(ee);
%     bLi = problem_data.bL(ee);
%     cLi = problem_data.cL(ee);
%     base_part = 2*cLi * lambda_e + bLi;
%     relu_part = 0.0;
%     if (ssom_relu_argument(lambda_e)>0)
%         relu_part = -rho;
%     end
%     g(ee) = base_part + relu_part;
% end
% end

function g=grad_lambda(lambdas, R, T, problem_data)

%%

edges = problem_data.edges;
tijs_vec = problem_data.tijs;
% rho = problem_data.rho;

num_edges = size(edges, 1);
aL = zeros(size(lambdas));
bL = zeros(size(lambdas));
cL = zeros(size(lambdas));

% cost_out = 0.0;
for ee = 1:num_edges
    ii = edges(ee, 1);
    jj = edges(ee, 2);
    % lambda_e = lambdas(ee);
    tij_e = tijs_vec(:, ee);
    T_i = T(:, ii);
    T_j = T(:, jj);
    R_i = R(:, :, ii);
    a = T_i - T_j;
    b = R_i * tij_e;
    % cost_lambda_0_ee = trace(aL' * aL + 2 * lambda_e * (aL' * b) + lambda_e^2 * (b' * b));
    % cost_relu_ee = relu_som(ssom_relu_argument(lambda_e));
    % cost_out = cost_out + cost_lambda_0_ee + rho * cost_relu_ee;
    aL(ee) = trace(a' * a);
    bL(ee) = 2 * trace(a' * b);
    cL(ee) = trace(b' * b);
end
%%


rho = problem_data.rho;

% x = lambdas in this context

g = zeros(length(lambdas), 1);

num_edges = size(edges, 1);
for ee = 1:num_edges
    % ii = edges(ee, 1);
    % jj = edges(ee, 2);
    lambda_e = lambdas(ee);
    % tij_e = Tijs_vec(:, ee);
    % aLi = problem_data.aL(ee);
    bLi = problem_data.bL(ee);
    cLi = problem_data.cL(ee);
    base_part = 2*cLi * lambda_e + bLi;
    relu_part = 0.0;
    if (ssom_relu_argument(lambda_e)>0)
        relu_part = -rho;
    end
    g(ee) = base_part + relu_part;
end
end

%% 1
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

%% 2
function h = ssom_ehess_r_lambda(R, Rdot, T, lambdadot, problem_data)

% h_r_lambda = zeros(size(hrt));

% lambdas = X.lambda;
% lambdas_dot = Xdot.lambda;
edges = problem_data.edges;
% Tijs_scaled = make_tijs_scaled(X.lambdas, problem_data.Tijs);
% rho = problem_data.rho;

W = zeros(size(R));
num_edges = size(edges, 1);
for e = 1:num_edges
    ii = edges(e,1);
    jj = edges(e,2);
    T_i = T(:,ii);
    T_j = T(:,jj);
    tij = problem_data.Tijs(:, e);
    lambda_dot_e = lambdadot(e);
    w_ij = 2 * (T_i - T_j)*lambda_dot_e*tij';
    W(:,:,ii) = W(:,:,ii) + w_ij;
end
h = stiefel_tangentProj(R, W);
% h = W;
end


%% 3
function h = ssom_ehess_t_lambda(x, ~, R, lambdas_dot, problem_data)
% h_t_lambda = zeros(size(htr));

edges = problem_data.edges;
% Tijs_scaled = make_tijs_scaled(X.lambdas, problem_data.Tijs);
% rho = problem_data.rho;

h = zeros(size(x));
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
    tij = problem_data.Tijs(:, e);
    w_ij = BIJ * lambdas_dot(e) * tij' * Ri';
    h = h + w_ij';
end

end



%% 4
function h = ssom_ehess_lambda_r(x, xdot, T, R, Rdot, problem_data)
% h_lambda_t = zeros(size(h_lambda_lambda));

% x = X.lambda;
% lambdas_dot = Xdot.lambda;
edges = problem_data.edges;
Tijs_vec = problem_data.Tijs;
% rho = problem_data.rho;

% nrs = problem_data.sz(1);
% % d = problem_data.sz(2);
% N = problem_data.sz(3);

h = zeros(size(x));

num_edges = size(edges, 1);
for ee = 1:num_edges
    ii = edges(ee, 1);
    jj = edges(ee, 2);
    lambda_e = x(ee);
    tij = Tijs_vec(:, ee);
    T_i = T(:, ii);
    T_j = T(:, jj);
    R_i_dot = Rdot(:, :, ii);
    R_i = R(:, :, ii);
    a = T_i - T_j;
    bdot = R_i_dot * tij;
    b = R_i * tij;
    e_th_elem = 2 * lambda_e * (bdot' * b + b' * bdot) + 2 * a' * bdot;
    h(ee) = e_th_elem;
end
end


%% 5
function h = ssom_ehess_lambda_t(x, ~, ~, Tdot, R, problem_data)

% h_lambda_r = zeros(size(h_lambda_lambda));

% x = X.lambda;
% lambdas_dot = Xdot.lambda;
edges = problem_data.edges;
tijs_vec = problem_data.tijs;
% rho = problem_data.rho;

h = zeros(size(x));

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


