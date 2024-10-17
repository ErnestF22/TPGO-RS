close all;
clear;
clc;


problem=test_problem();
curve=test_problem_curve(problem);


manif = stiefelfactory(problem.sz(1),problem.sz(2),problem.sz(3));
problem_manopt.M = manif;

problem_manopt.cost = @(x) som_cost_rot_stiefel(x, problem);
problem_manopt.egrad = @(x) som_egrad_rot_stiefel(x, problem);
problem_manopt.grad = @(x) som_rgrad_rot_stiefel(x, problem);
problem_manopt.ehess = @(x,u) som_ehess_rot_stiefel(x,u, problem);
problem_manopt.hess = @(x,u) som_rhess_rot_stiefel(x,u, problem);


% R_initguess = eye3d(problem.sz(1),problem.sz(2),problem.sz(3));
R_initguess = make_rand_stiefel_3d_array(problem.sz(1),problem.sz(2),problem.sz(3));
options.maxiter = 100;
[R, R_cost, R_info, R_options] = trustregions(problem_manopt, R_initguess, options);

disp("problem_manopt.grad(R)")
disp(problem_manopt.grad(R))
disp("max(problem_manopt.grad(R))");
disp(max(problem_manopt.grad(R),[], 'all'));
L = problem.L;
P = problem.P;
[lambda_pim, v_pim] = pim_hessian(R, problem);

R_next = cat_zero_rows_3d_array(R);
% v_pim_next = cat_zero_rows_3d_array(v_pim);
problem_next = problem;
problem_next.sz(1) = problem_next.sz(1) + 1;
L_next = from_L_to_L_next(L,problem_next);
P_next = from_P_to_P_next(P,problem_next);
problem_next.L = L_next;
problem_next.P = P_next;
[lambda_pim_next, v_pim_next] = pim_hessian(R_next, problem_next);

%Plot

%!! the function plotted (after padding with zeros) will need to be concave
%linesearch will give the new point with decreased cost (i.e., the initial
%guess point for the next step of the Riemannian Staircase)

lambdas = reshape(-1:0.1:1, 1,1,size(-1:0.1:1, 2));
lambdas_to_add = lambdas + repmat(matStack(v_pim), 1,1, size(lambdas,3));
x_plot = repmat(matStack(R), 1,1, size(lambdas,3)) + lambdas_to_add;
fx_plot = zeros(size(lambdas, 3), 1);
for ii = 1:size(lambdas, 3)
    matunst_x = matUnstack(x_plot(:,:,ii), size(R,1));
    fx_plot(ii) = som_cost_rot_stiefel(matunst_x, problem);
end

% x_plot_plot
x_plot_plot = zeros(problem.sz(1)*problem.sz(2)*problem.sz(3), size(lambdas,3));
for ii = 1:size(lambdas, 3)
    x_plot_ii = x_plot(:,:,ii);
    x_plot_plot(:,ii) = x_plot_ii(:);
end
% x_plot_plot = x_plot_plot(1,:);
plot(x_plot_plot, fx_plot)


