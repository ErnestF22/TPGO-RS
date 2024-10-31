function scale_formulation_manopt

% importmanopt;


nrs = 4;
d = 3;
N = 5;



num_edges = 8;
G = graph(true(N), 'omitselfloops'); % Alternative without self-loops
p = randperm(numedges(G), num_edges);
G = graph(G.Edges(p, :));
edges = table2array(G.Edges);


Tijs_vec = 10 * rand(d, num_edges);
problem_data.Tijs = Tijs_vec;
problem_data.edges = edges;


tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
tuple.lambda = euclideanfactory(num_edges, 1);
M = productmanifold(tuple);

problem_data.rho = 0.0;

problem.M = M;
problem_data.sz = [nrs, d, N];
problem.cost = @(x) ssom_cost(x, problem_data);
problem.egrad = @(x) ssom_egrad(x, problem_data);
% problem.grad = @(x) ssom_rgrad(x, problem_data);
problem.hess = @(x, u) ssom_hess_genproc(x, u, problem_data);

fprintf("\n");
disp("Gradient check easy");
X_gradcheck.R = eye3d(nrs, d, N);
X_gradcheck.T = zeros(nrs, N);
X_gradcheck.lambda = ones(num_edges, 1);
figure(1)
checkgradient(problem, X_gradcheck);

fprintf("\n");
disp("Gradient check rand");
figure(2)
checkgradient(problem);

close all;

fprintf("\n");
disp("Hessian check");
figure(3)
checkhessian(problem);

% problem_data.R_gt = eye3d(nrs, d, N);
% problem_data.T_gt = zeros(nrs, N);
% problem_data.lambda_gt = zeros(num_edges, 1);



% trustregions(problem);


end %file function

%%
function [h] = ssom_hess_genproc(X, Xdot, problem_data)
    R = X.R;
    T = X.T;
    lambdas = X.lambda;
    Rdot = Xdot.R;
    Tdot = Xdot.T;
    lambdas_dot = Xdot.lambda;
    
    % Careful: tangent vectors on the rotation group are represented as
    % skew symmetric matrices. To obtain the corresponding vectors in
    % the ambient space, we need a little transformation. This
    % transformation is typically not needed when we compute the
    % formulas for the gradient and the Hessian directly in Riemannian
    % form instead of resorting the egrad2rgrad and ehess2rhess. These
    % latter tools are convenient for prototyping but are not always
    % the most efficient form to execute the computations.
    edges = problem_data.edges;

    tijs_scaled = make_tijs_scaled(lambdas, problem_data.Tijs);
    
    %g.R
    problem_data_R = problem_data;
    problem_data_R.Tijs = tijs_scaled;
    [problem_data_R.P, problem_data_R.frct] = ...
        make_step1_p_fct(T, tijs_scaled, edges);
    g.R = rsom_rgrad_rot_stiefel(R, problem_data_R);
    %g.T
    problem_data_T = problem_data;
    problem_data_T.Tijs = tijs_scaled;
    [problem_data_T.LR, problem_data_T.PR, problem_data_T.BR] = ...
        make_LR_PR_BR_noloops(R, tijs_scaled, edges);
    g.T = rsom_rgrad_transl_stiefel(T, problem_data_T);
    %g.lambda
    problem_data_lambdas = problem_data;
    problem_data_lambdas.T = T;
    problem_data_lambdas.R = R;
    g.lambda = ssom_grad_lambda(lambdas, problem_data_lambdas);
 
    hrt = compute_hrt(X,Xdot,tijs_scaled,edges);
    htr = compute_htr(X,Xdot,tijs_scaled,edges);
    
    % h_lambda_lambda = zeros(1,1); %TODO1
    h_lambda_lambda = ssom_rhess_lambda(X, Xdot, problem_data_lambdas);

    h_r_lambda = zeros(size(hrt)); %TODO 2
    h_t_lambda = zeros(size(htr)); %TODO 3

    h_lambda_r = zeros(size(h_lambda_lambda)); %TODO 4
    h_lambda_t = zeros(size(h_lambda_lambda)); %TODO 5

    h.R = rsom_rhess_rot_stiefel(R, Rdot, problem_data_R) + hrt + h_r_lambda;
    h.T = rsom_rhess_transl_stiefel(T, Tdot, problem_data_T) + htr + h_t_lambda;
    h.lambda = h_lambda_lambda + h_lambda_r + h_lambda_t;
end %hess

%%
function h = ssom_rhess_lambda(X, Xdot, problem_data)
    lambdas = X.lambda;
    lambdas_dot = Xdot.lambda;
    edges = problem_data.edges;
    Tijs_vec = problem_data.Tijs;
    % rho = problem_data.rho;


    h = zeros(length(lambdas), 1);
    
    num_edges = size(edges, 1);
    for ee = 1:num_edges
        ii = edges(ee, 1);
        % jj = edges(ee, 2);
        % lambda_e = lambdas(ee);
        tij_e = Tijs_vec(:, ee);
        % T_i = problem_data.T(:, ii);
        % T_j = problem_data.T(:, jj);
        R_i = problem_data.R(:, :, ii);
        % a = T_i - T_j;
        b = R_i * tij_e;
        base_part = 2*(b' * b);        
        h(ee) = base_part * lambdas_dot(ee);
    end
    
end