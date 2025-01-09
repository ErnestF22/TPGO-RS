function scale_formulation_manopt

% importmanopt;


nrs = 4;
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


rho = 1.0; %TODO: make this rand() later

% variables random generation/init
tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
tuple.lambda = euclideanfactory(num_edges, 1);
M = productmanifold(tuple);
% problem.M = M;
% X = M.rand();

problem_data = struct('sz', sz, 'edges', edges, 'tijs', tijs);

%g.lambda
% [aL, bL, cL] = makeABClambda(X, problem_data);
% problem_data.aL = aL; problem_data.bL = bL; problem_data.cL = cL;

problem_data.rho = rho; %ReLU() part should not be needed for Hessian tests

%% output (problem_curve_data) definition


problem.M = M;
problem_data.sz = [nrs, d, N];
problem.cost = @(x) ssom_cost(x, problem_data);
problem.egrad = @(x) ssom_egrad(x, problem_data);
% problem.grad = @(x) ssom_rgrad(x, problem_data);
problem.ehess = @(x, u) ssom_ehess_genproc(x, u, problem_data);

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
disp("Hessian check");
% X_gradcheck.R = eye3d(nrs, d, N);
% X_gradcheck.T = zeros(nrs, N);
% X_gradcheck.lambda = zeros(num_edges, 1);
% Xdot_gradcheck.R = eye3d(nrs, d, N);
% Xdot_gradcheck.T = zeros(nrs, N);
% Xdot_gradcheck.lambda = zeros(num_edges, 1);
figure(3)
checkhessian(problem);

% problem_data.R_gt = eye3d(nrs, d, N);
% problem_data.T_gt = zeros(nrs, N);
% problem_data.lambda_gt = zeros(num_edges, 1);



% trustregions(problem);


end %file function

%%
function [h] = ssom_ehess_genproc(X, Xdot, problem_data)
    R = X.R;
    T = X.T;
    lambdas = X.lambda;
    Rdot = Xdot.R;
    Tdot = Xdot.T;
    % lambdas_dot = Xdot.lambda;
    
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
    
    %g.R
    problem_data_R = problem_data;
    problem_data_R.Tijs = tijs_scaled;
    [problem_data_R.P, problem_data_R.frct] = ...
        make_step1_p_fct(T, tijs_scaled, edges);
    % g.R = rsom_rgrad_rot_stiefel(R, problem_data_R);
    %g.T
    problem_data_T = problem_data;
    problem_data_T.Tijs = tijs_scaled;
    [problem_data_T.LR, problem_data_T.PR, problem_data_T.BR] = ...
        make_LR_PR_BR_noloops(R, tijs_scaled, edges);
    % g.T = rsom_rgrad_transl_stiefel(T, problem_data_T);
    %g.lambda
    problem_data_lambdas = problem_data;
    problem_data_lambdas.T = T;
    problem_data_lambdas.R = R;
    % g.lambda = ssom_grad_lambda(lambdas, problem_data_lambdas);
 
    hrt = compute_hrt(X,Xdot,tijs_scaled,edges);
    htr = compute_htr(X,Xdot,tijs_scaled,edges);
    
    % h_lambda_lambda = zeros(1,1); %TODO1
    h_lambda_lambda = ssom_ehess_lambda_lambda(X, Xdot, problem_data_lambdas);

    % h_r_lambda = zeros(size(hrt)); %TODO 2
    h_r_lambda = ssom_ehess_r_lambda(X, Xdot, problem_data_lambdas);

    % h_t_lambda = zeros(size(htr)); %TODO 3
    h_t_lambda = ssom_ehess_t_lambda(X, Xdot, problem_data_lambdas);

    % h_lambda_r = zeros(size(h_lambda_lambda)); %TODO 4
    h_lambda_r = ssom_ehess_lambda_r(X, Xdot, problem_data_lambdas);
    
    % h_lambda_t = zeros(size(h_lambda_lambda)); %TODO 5
    h_lambda_t = ssom_ehess_lambda_t(X, Xdot, problem_data_lambdas);

    h.R = rsom_ehess_rot_stiefel(R, Rdot, problem_data_R) + hrt + h_r_lambda;
    h.T = rsom_ehess_transl_stiefel(T, Tdot, problem_data_T) + htr + h_t_lambda;
    h.lambda = h_lambda_lambda + h_lambda_r + h_lambda_t;
end %hess

%% 1
function h = ssom_ehess_lambda_lambda(X, Xdot, problem_data)
    % h_lambda_lambda = zeros(1,1)

    lambdas = X.lambda;
    % lambdas_dot = Xdot.lambda;
    edges = problem_data.edges;
    tijs = problem_data.tijs;
    % rho = problem_data.rho;

    h = zeros(length(lambdas), 1);
    
    num_edges = size(edges, 1);
    for ee = 1:num_edges
        ii = edges(ee, 1);
        % jj = edges(ee, 2);
        % lambda_e = lambdas(ee);
        tij_e = tijs(:, ee);
        % T_i = problem_data.T(:, ii);
        % T_j = problem_data.T(:, jj);
        R_i = X.R(:, :, ii);
        % a = T_i - T_j;
        b = R_i * tij_e;    
        h(ee) = 2*(b' * b) * Xdot.lambda(ee);
    end
    
end

%% 2
function h = ssom_ehess_r_lambda(X, Xdot, problem_data)
    % h_r_lambda = zeros(size(hrt));

    % lambdas = X.lambda;
    lambdas_dot = Xdot.lambda;
    edges = problem_data.edges;
    % rho = problem_data.rho;

    tijs_scaled = make_tijs_scaled(X.lambda, problem_data.tijs);

    W = zeros(size(X.R));
    num_edges = size(edges, 1);
    for e = 1:num_edges
        ii = edges(e,1);
        jj = edges(e,2); 
        T_i = X.T(:,ii);
        T_j = X.T(:,jj);
        tij = tijs_scaled(:, e);
        lambda_dot_e = lambdas_dot(e);
        w_ij = 2 * (T_i - T_j)*lambda_dot_e*tij';
        W(:,:,ii) = W(:,:,ii) + w_ij;
    end
    h = stiefel_tangentProj(X.R, W);
    % h = W;
end

%% 3
function h = ssom_ehess_t_lambda(X, Xdot, problem_data)
    % h_t_lambda = zeros(size(htr));

    lambdas_dot = Xdot.lambda;
    edges = problem_data.edges;
    tijs_scaled = make_tijs_scaled(X.lambda, problem_data.tijs);
    % rho = problem_data.rho;

    h = zeros(size(X.T));
    N = size(X.R,3);
    num_edges = size(edges, 1);
    for e = 1:num_edges
        ii = edges(e,1);
        jj = edges(e,2); 
        Ri = X.R(:,:,ii);
%         Rj = X.R(:,:,jj);
        %
        BIJ = zeros(N,1);
        BIJ(ii) = 1;
        BIJ(jj) = -1;
        %
        tij = tijs_scaled(:, e);
        w_ij = BIJ * lambdas_dot(e) * tij' * Ri';
        h = h + w_ij';
    end

end

%% 4
function h = ssom_ehess_lambda_r(X, Xdot, problem_data)
    % h_lambda_r = zeros(size(h_lambda_lambda));

    lambdas = X.lambda;
    % lambdas_dot = Xdot.lambda;
    edges = problem_data.edges;
    tijs = problem_data.tijs;
    % tijs_scaled = make_tijs_scaled(X.lambdas, problem_data.tijs);
    % rho = problem_data.rho;

    h = zeros(size(lambdas));
    
    num_edges = size(edges, 1);
    for ee = 1:num_edges
        ii = edges(ee, 1);
        jj = edges(ee, 2);
        lambda_e = lambdas(ee);
        tij = tijs(:, ee);
        T_i = X.T(:,ii);
        T_j = X.T(:,jj);
        T_i_dot = Xdot.T(:, ii);
        T_j_dot = Xdot.T(:, jj);
        a = T_i - T_j;
        R_i = X.R(:, :, ii);
        b = R_i * tij;
        adot = T_i_dot - T_j_dot;
        e_th_elem = adot' * a + a' * adot + 2 * lambda_e * adot' * b;
        h(ee) = e_th_elem;
    end

end

%% 5
function h = ssom_ehess_lambda_t(X, Xdot, problem_data)
    % h_lambda_t = zeros(size(h_lambda_lambda));

    lambdas = X.lambda;
    % lambdas_dot = Xdot.lambda;
    edges = problem_data.edges;
    Tijs_vec = problem_data.tijs;
    % rho = problem_data.rho;

    % nrs = problem_data.sz(1);
    % % d = problem_data.sz(2);
    % N = problem_data.sz(3);
    
    h = zeros(size(lambdas));
    
    num_edges = size(edges, 1);
    for ee = 1:num_edges
        ii = edges(ee, 1);
        jj = edges(ee, 2);
        lambda_e = lambdas(ee);
        tij = Tijs_vec(:, ee);
        T_i = X.T(:, ii);
        T_j = X.T(:, jj);
        R_i_dot = Xdot.R(:, :, ii);
        R_i = X.R(:, :, ii);
        a = T_i - T_j;
        bdot = R_i_dot * tij;
        b = R_i * tij;
        e_th_elem = lambda_e^2 * (bdot' * b + b' * bdot) + 2 * lambda_e * a' * bdot;
        h(ee) = e_th_elem;
    end

end
%%
