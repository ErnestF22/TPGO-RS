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

%%%%%
M1 = make_rand_stiefel_3d_array(problem_struct.sz(1),problem_struct.sz(2),problem_struct.sz(3))
M2 = make_rand_stiefel_3d_array(problem_struct.sz(1),problem_struct.sz(2),problem_struct.sz(3))

problem_manopt.M.norm(M1, R)
problem_manopt.M.norm([], R)
problem_manopt.M.norm([], M2)
problem_manopt.M.norm([], [])
problem_manopt.M.norm(M2, M1)
norm(R(:))
norm(M2(:))

%%%%%
disp()
problem_manopt.M.inner(M1, R,R)
problem_manopt.M.inner([], R, R)
problem_manopt.M.inner([], M2, M2)
problem_manopt.M.inner([], [], [])
problem_manopt.M.inner(M2, M1, R)


