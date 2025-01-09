function h = ssom_rhess_genproc(X, Xdot, problem_data)
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
    edges = problem_data.edges;

    tijs_scaled = make_tijs_scaled(lambdas, problem_data.tijs);
    
    % problem helpful members
    [problem_data.P, problem_data.frct] = ...
        make_step1_p_fct(T, tijs_scaled, edges);
    %
    [problem_data.LR, problem_data.PR, problem_data.BR] = ...
        make_LR_PR_BR_noloops(R, tijs_scaled, edges);
    %
 
    hrt = compute_hrt(X,Xdot,tijs_scaled,edges);
    htr = compute_htr(X,Xdot,tijs_scaled,edges);
    
    % h_lambda_lambda = zeros(size(lambda));
    h_lambda_lambda = ssom_ehess_lambda_lambda(lambdas, lambdasdot, R, problem_data);

    % h_r_lambda = zeros(size(hrt)); 
    h_r_lambda = ssom_ehess_r_lambda(R, T, lambdas, lambdasdot, problem_data);

    % h_t_lambda = zeros(size(htr)); 
    h_t_lambda = ssom_ehess_t_lambda(R, T, lambdas, lambdasdot, problem_data);

    % h_lambda_r = zeros(size(h_lambda_lambda)); 
    h_lambda_r = ssom_rhess_lambda_r(X, Xdot, problem_data);
    
    % h_lambda_t = zeros(size(h_lambda_lambda)); 
    h_lambda_t = ssom_ehess_lambda_t(R, T, Tdot, lambdas, problem_data);

    h.R = ssom_rhess_R_R(R, Rdot, problem_data) + hrt + h_r_lambda;
    h.T = ssom_ehess_T_T(T, Tdot, problem_data) + htr + h_t_lambda;
    h.lambda = h_lambda_lambda + h_lambda_r + h_lambda_t;
    
end

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

function h = ssom_ehess_T_T(~, xdot, problem_data)
h = xdot*(problem_data.LR' + problem_data.LR);
end

function h = ssom_rhess_R_R(x, xdot, problem_data)
d = size(x, 2);
egrad = matUnstackH(problem_data.P,d); %!! ehess2rhess for stiefel manifolds!
h = ehess2rhess_stiefel(x, xdot, egrad);
end

function rhess = ehess2rhess_stiefel(x, xdot, egrad)
term_1 = multiprod(xdot, ...
    0.5*multiprod(multitransp(x), egrad) + 0.5*multiprod(multitransp(egrad), x));
term_2 = multiprod(x, ...
    0.5*multiprod(multitransp(xdot), egrad) + 0.5*multiprod(multitransp(egrad), xdot));
DGf = - term_1 - term_2;
rhess = stiefel_tangentProj(x, DGf); %ehess_proj = zeros(nrs,d,N)
end

%% 2
function h = ssom_ehess_r_lambda(R, T, ~, lambdadot, problem_data)

% h_r_lambda = zeros(size(hrt));

% lambdas = X.lambda;
% lambdas_dot = Xdot.lambda;
edges = problem_data.edges;
% rho = problem_data.rho;

W = zeros(size(R));
num_edges = size(edges, 1);
for e = 1:num_edges
    ii = edges(e,1);
    jj = edges(e,2);
    T_i = T(:,ii);
    T_j = T(:,jj);
    tij = problem_data.tijs(:, e);
    lambda_dot_e = lambdadot(e);
    w_ij = 2 * (T_i - T_j)*lambda_dot_e*tij';
    W(:,:,ii) = W(:,:,ii) + w_ij;
end
h = stiefel_tangentProj(R, W);
% h = W;
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
    h = h + w_ij';
end



end



%% 4 
function eh = ssom_ehess_lambda_r(R, Rdot, T, lambdas, problem_data)
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
    e_th_elem = 2 * lambda_e * (bdot' * b + b' * bdot) + 2 * a' * bdot;
    eh(ee) = e_th_elem;
end

end

function h = ssom_rhess_lambda_r(X, Xdot, problem_data)
R = X.R;
Rdot = Xdot.R;
T = X.T;
lambdas = X.lambda;
lambdasdot = Xdot.lambda;
eh = ssom_ehess_lambda_r(R,Rdot,T,lambdas, problem_data);

%HP) H = u

eg = ssom_grad_lambda(X, problem_data);
h = manopt_stiefel_ehess2rhess(lambdas, eg, eh, lambdasdot);

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

function rhess = manopt_stiefel_ehess2rhess(X, egrad, ehess, H)
    XtG = multiprod(multitransp(X), egrad);
    symXtG = multisym(XtG);
    HsymXtG = multiprod(H, symXtG);
    rhess = stiefel_tangentProj(X, ehess - HsymXtG);
end
