function test_problem_test_pim_nextstep

problem_struct = test_problem();
curve=test_problem_curve(problem_struct);


manif = stiefelfactory(problem_struct.sz(1),problem_struct.sz(2),problem_struct.sz(3));
problem_manopt.M = manif;

problem_manopt.cost = @(x) som_cost_rot_stiefel(x, problem_struct);
problem_manopt.egrad = @(x) som_egrad_rot_stiefel(x, problem_struct);
problem_manopt.rgrad = @(x) som_rgrad_rot_stiefel(x, problem_struct);
problem_manopt.ehess = @(x,u) som_ehess_rot_stiefel(x,u, problem_struct);
problem_manopt.rhess = @(x,u) som_rhess_rot_stiefel(x,u, problem_struct);


R_initguess = eye3d(problem_struct.sz(1),problem_struct.sz(2),problem_struct.sz(3));
options.maxiter = 100;
[R, R_cost, R_info, R_options] = trustregions(problem_manopt, R_initguess, options);

fprintf("\nproblem_manopt.rgrad(R)\n")
disp(problem_manopt.rgrad(R))
fprintf("\nmax(problem_manopt.rgrad(R))\n")
disp(max(problem_manopt.rgrad(R),[], 'all'));

%Run P.I.M. on Manopt's output
[lambda_pim, v_pim] = pim_hessian(R, problem_struct);

%check if lambda_pim is actually an eigenvalue:
disp("is lambda_pim an eigenvalue?")
disp([lambda_pim * v_pim, som_rhess_rot_stiefel(R, v_pim, problem_struct)])

R_next = cat_zero_rows_3d_array(R);
v_pim_next = cat_zero_rows_3d_array(v_pim);
problem_next_struct = problem_struct;
problem_next_struct.sz(1) = problem_next_struct.sz(1) + 1;
L_next = from_L_to_L_next(problem_struct.L, problem_next_struct);
P_next = from_P_to_P_next(problem_struct.P, problem_next_struct);
problem_next_struct.L = L_next;
problem_next_struct.P = P_next;
[lambda_pim_next, v_pim_next] = pim_hessian(R_next, problem_next_struct);



disp("is lambda_pim_next an eigenvalue?")
disp([lambda_pim_next * v_pim_next, ...
    som_rhess_rot_stiefel(R_next, v_pim_next, problem_next_struct)])

end %file function
