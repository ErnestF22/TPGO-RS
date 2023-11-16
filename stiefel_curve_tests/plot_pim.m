problem=test_problem();
curve=test_problem_curve(problem);


manif = stiefelfactory(problem.sz(1),problem.sz(2),problem.sz(3));
problem_manopt.M = manif;

problem_manopt.cost = @(x) som_cost_rot_stiefel(x, problem);
problem_manopt.egrad = @(x) som_egrad_rot_stiefel(x, problem);
problem_manopt.rgrad = @(x) som_rgrad_rot_stiefel(x, problem);
problem_manopt.ehess = @(x,u) som_ehess_rot_stiefel(x,u, problem);
problem_manopt.rhess = @(x,u) som_rhess_rot_stiefel(x,u, problem);


R_initguess = eye3d(problem.sz(1),problem.sz(2),problem.sz(3));
options.maxiter = 100;
[R, R_cost, R_info, R_options] = trustregions(problem_manopt, R_initguess, options);

problem_manopt.rgrad(R)
max(problem_manopt.rgrad(R),[], 'all')
L = @(u) problem_manopt.rhess(R,u);
[lambda_pim, v_pim] = pim_hessian(R, problem);

R_next = cat_zero_rows_3d_array(R);
% v_pim_next = cat_zero_rows_3d_array(v_pim);
problem_next = problem;
L_next = @(u) problem_manopt.rhess(R_next,u);
problem_next.sz(1) = problem_next.sz(1) + 1;
[lambda_pim_next, v_pim_next] = pim_hessian(R_next, L_next, problem_next);

%plot...

%!! the function plotted (after padding with zeros) will need to be concave
%linesearch will give the new point with decreased cost (i.e., the initial
%guess point for the next step of the Riemannian Staircase)


