function [X_manopt_out] = lsom_rtr(nrs, d, N, problem_data, params, transf_initguess_struct, lambdas_initguess)

edges = problem_data.E;

num_edges = size(edges, 1);

tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
tuple.lambda = euclideanfactory(num_edges, 1);
M = productmanifold(tuple);

% Setup the problem structure with manifold M and cost+grad functions.
problem.M = M;

if params.relu_scale_compensation
    problem.cost = @(x) lsom_cost_relu(x, problem_data);
    % problem.egrad = @(x) lsom_egrad(x, problem_data);
    problem.grad = @(x) lsom_rgrad_relu(x, problem_data);
    % problem.ehess = @(x, u) lsom_ehess_genproc(x, u, problem_data);
    problem.hess = @(x, u) lsom_rhess_genproc_relu(x, u, problem_data);
else
    problem.cost = @(x) lsom_cost(x, problem_data);
    % problem.egrad = @(x) lsom_egrad(x, problem_data);
    problem.grad = @(x) lsom_rgrad(x, problem_data);
    % problem.ehess = @(x, u) lsom_ehess_genproc(x, u, problem_data);
    problem.hess = @(x, u) lsom_rhess_genproc(x, u, problem_data);
end

% checkgradient(problem);
% tmp.R = make_rand_stiefel_3d_array(nrs, d, N);
% tmp.R = eye3d(nrs, d, N);
% tmp.T = rand(nrs, N);
% tmp.lambda = rand(num_edges, 1);
% tmpU.R = zeros(nrs, d, N);
% tmpU.T = rand(nrs, N);
% tmpU.T = normalize(tmpU.T);
% tmpU.lambda = rand(num_edges, 1);
% tmpU.lambda = normalize(tmpU.lambda);
close all;
figure(10)
% checkgradient(problem, tmp);
X_chkgrad = M.rand();
X_chkgrad.lambda = 2 * problem_data.a + rand(num_edges, 1); %!! function not defined on all lambdas
X_chkgrad_tg = M.randvec(X_chkgrad);
% X_chkgrad_tg.lambda = 2 * problem_data.a + rand(num_edges, 1);
checkgradient(problem, X_chkgrad,X_chkgrad_tg)
figure(11)
% checkhessian(problem, tmp);
checkhessian(problem, X_chkgrad, X_chkgrad_tg)

%check that GT cost is 0
% !! only works when tijs are gt
X_gt.lambda = problem_data.lambda_gt;
X_gt.R = problem_data.R_gt;
X_gt.T = problem_data.T_gt;

if params.relu_scale_compensation
    cost_gt = lsom_cost_relu(X_gt, problem_data);
    disp("cost_gt_relu in lsom_genproc.m")
    disp(cost_gt)
else
    cost_gt = lsom_cost(X_gt, problem_data);
    disp("cost_gt in lsom_genproc.m")
    disp(cost_gt)
end

disp("cost gt _no_compensation(X_recovered, problem_data)")
disp(ssom_cost_no_compensation(X_gt, problem_data))

% tg_element_test = M.randvec(X_gt);
% disp("check_is_tangent_stiefel(X_gt.R, tg_element_test.R)")
% disp(check_is_tangent_stiefel(X_gt.R, tg_element_test.R));


% X = trustregions(problem, X_gt);
options.maxiter = 1000;

X_initguess.R = transf_initguess_struct.R;
X_initguess.T = transf_initguess_struct.T;
X_initguess.lambda = lambdas_initguess;

% disp("transf_initguess")
% disp(transf_initguess)


if params.relu_scale_compensation
    cost_initguess = lsom_cost_relu(X_initguess, problem_data);
    disp("cost_initguess")
    disp(cost_initguess)
else
    cost_initguess = lsom_cost(X_initguess, problem_data);
    disp("cost_initguess")
    disp(cost_initguess)
end

% rg_ig = lsom_rgrad(X_initguess, problem_data);
% disp("lsom_rgrad(X_initguess, problem_data)")
% disp(rg_ig.R)
% disp(rg_ig.T)
% disp(rg_ig.lambda)

X = trustregions(problem, X_initguess, options);
T_manopt_out = X.T;
R_manopt_out = X.R;
lambdas_manopt_out = X.lambda;




X_manopt_out.R = R_manopt_out;
X_manopt_out.T = T_manopt_out;
X_manopt_out.lambda = lambdas_manopt_out;

if params.relu_scale_compensation
    cost_manopt_out = lsom_cost_relu(X_manopt_out, problem_data);
else
    cost_manopt_out = lsom_cost(X_manopt_out, problem_data);
end
disp("cost_manopt_out")
disp(cost_manopt_out)

end %file function