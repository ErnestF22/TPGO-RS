function problem_curve_data=test_check_ssom()

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

Tijs = 10 * rand(d, num_edges);

problem_data = struct('sz', sz, 'edges', edges, 'Tijs', Tijs);

lambdas = 10 * rand(num_edges, 1);

rho = 1.0; %TODO: make this rand() later

tijs_scaled = make_tijs_scaled(lambdas, problem_data.Tijs); %!!

% variables random generation/init
tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
tuple.lambda = euclideanfactory(num_edges, 1);
M = productmanifold(tuple);
problem.M = M;
problem_data.sz = [nrs, d, N];
x = M.rand();
% lambdas = x.lambda;
T = x.T;
R = x.R;

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
[aL, bL, cL] = makeABClambda(x, problem_data_lambdas);



problem.sz = sz;
problem.P = problem_data_R.P; problem.frct = problem_data_R.frct;
problem.LR = problem_data_T.LR; problem.PR = problem_data_T.PR; problem.BR = problem_data_T.BR;
problem.rho = rho; %ReLU() part should not be needed for Hessian tests
problem.aL = aL; problem.bL = bL; problem.cL = cL;

%% output (problem_curve_data) definition

problem_curve_data = problem;
problem_curve_data.edges = edges;
problem_curve_data.Tijs = Tijs;
problem_curve_data.cost=@(x) ssom_cost_lambda(x,problem_curve_data);
problem_curve_data.egrad_lambda=@(x) egrad_lambda(x,problem_curve_data);
problem_curve_data.ssom_ehess_lambda_lambda=@(x,u) ssom_ehess_lambda_lambda(x,u,problem_curve_data);
problem_curve_data.ssom_ehess_lambda_r=@(x,u) ssom_ehess_lambda_r(x,u,problem_curve_data);
problem_curve_data.ssom_ehess_r_lambda=@(x,u) ssom_ehess_r_lambda(x,u,problem_curve_data);
problem_curve_data.ssom_ehess_lambda_t=@(x,u) ssom_ehess_lambda_t(x,u,problem_curve_data);
problem_curve_data.ssom_ehess_t_lambda=@(x,u) ssom_ehess_t_lambda(x,u,problem_curve_data);

end %file function


%% 0) egrad lambda
function g=egrad_lambda(x,problem_data_lambdas)
    %Note: already checked (genproc-wise) through Manopt
    
    % x in this context is already lambda

    % g = ssom_grad_lambda(x, problem_data_lambdas);

    edges = problem_data_lambdas.edges;
    % Tijs_vec = problem_data_lambdas.Tijs;
    rho = problem_data_lambdas.rho;

    % x = lambdas in this context

    g = zeros(length(x), 1);
    
    num_edges = size(edges, 1);
    for ee = 1:num_edges
        % ii = edges(ee, 1);
        % jj = edges(ee, 2);
        lambda_e = x(ee);
        % tij_e = Tijs_vec(:, ee);
        % aLi = problem_data_lambdas.aL(ee);
        bLi = problem_data_lambdas.bL(ee);
        cLi = problem_data_lambdas.cL(ee);
        base_part = cLi * lambda_e + bLi;
        relu_part = 0.0;
        if (ssom_relu_argument(lambda_e)>0)
            relu_part = ssom_relu_argument(lambda_e);
        end
        g(ee) = base_part + rho * relu_part;
    end
end

%% 1
function h = ssom_ehess_lambda_lambda(X, Xdot, problem_data)
    % h_lambda_lambda = zeros(1,1)

    lambdas = X.lambda;
    % lambdas_dot = Xdot.lambda;
    edges = problem_data.edges;
    Tijs = problem_data.Tijs;
    % rho = problem_data.rho;

    h = zeros(length(lambdas), 1);
    
    num_edges = size(edges, 1);
    for ee = 1:num_edges
        ii = edges(ee, 1);
        % jj = edges(ee, 2);
        % lambda_e = lambdas(ee);
        tij_e = Tijs(:, ee);
        % T_i = problem_data.T(:, ii);
        % T_j = problem_data.T(:, jj);
        R_i = X.R(:, :, ii);
        % a = T_i - T_j;
        b = R_i * tij_e;    
        h(ee) = 2*(b' * b) * Xdot.lambda(ee);
    end
    
end

%% 2
function h = ssom_ehess_lambda_r(X, Xdot, problem_data)
    % h_lambda_r = zeros(size(h_lambda_lambda));

    lambdas = X.lambda;
    % lambdas_dot = Xdot.lambda;
    edges = problem_data.edges;
    Tijs_vec = problem_data.Tijs;
    % rho = problem_data.rho;

    h = zeros(size(lambdas));
    
    num_edges = size(edges, 1);
    for ee = 1:num_edges
        ii = edges(ee, 1);
        jj = edges(ee, 2);
        lambda_e = lambdas(ee);
        tij = Tijs_vec(:, ee);
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

%% 3
function h = ssom_ehess_r_lambda(X, Xdot, problem_data)
    % h_r_lambda = zeros(size(hrt));

    % lambdas = X.lambda;
    lambdas_dot = Xdot.lambda;
    edges = problem_data.edges;
    % Tijs_scaled = make_tijs_scaled(X.lambdas, problem_data.Tijs);
    % rho = problem_data.rho;

    W = zeros(size(X.R));
    num_edges = size(edges, 1);
    for e = 1:num_edges
        ii = edges(e,1);
        jj = edges(e,2); 
        T_i = X.T(:,ii);
        T_j = X.T(:,jj);
        tij = problem_data.Tijs(:, e);
        lambda_dot_e = lambdas_dot(e);
        w_ij = 2 * (T_i - T_j)*lambda_dot_e*tij';
        W(:,:,ii) = W(:,:,ii) + w_ij;
    end
    h = stiefel_tangentProj(X.R, W);
    % h = W;
end

%% 4
function h = ssom_ehess_lambda_t(X, Xdot, problem_data)
    % h_lambda_t = zeros(size(h_lambda_lambda));

    lambdas = X.lambda;
    % lambdas_dot = Xdot.lambda;
    edges = problem_data.edges;
    Tijs_vec = problem_data.Tijs;
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

%% 5
function h = ssom_ehess_t_lambda(X, Xdot, problem_data)
    % h_t_lambda = zeros(size(htr));

    lambdas_dot = Xdot.lambda;
    edges = problem_data.edges;
    % Tijs_scaled = make_tijs_scaled(X.lambdas, problem_data.Tijs);
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
        tij = problem_data.Tijs(:, e);
        w_ij = BIJ * lambdas_dot(e) * tij' * Ri';
        h = h + w_ij';
    end

end
