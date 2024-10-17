%% test call
test_rsom_pim_hessian


%% RSOM_PIM_HESSIAN  function call
[Y0_fun, lambda_pim_out_fun, v_pim_out_fun] = rsom_pim_hessian(R, problem_struct_next, thresh);
disp("[Y0, Y0_fun]");
disp([Y0, Y0_fun]);
disp("[lambda_pim_out_fun, lambda_pim_out_fun]");
disp([lambda_pim_out_fun, lambda_pim_out_fun]);
disp("[v_pim_out_fun, v_pim_out_fun]");
disp([v_pim_out_fun, v_pim_out_fun]);


%%
%%% check cost progression

disp("first RTR output (R) -> cost")
c_R = rsom_cost_rot_stiefel(R, problem_struct);
disp(c_R)

disp("first RTR output padded (x) -> cost")
c_x = rsom_cost_rot_stiefel(x, problem_struct_next);
disp(c_x)

disp("second RTR input (Y0) -> cost (FUN)")
c_Y0_fun = rsom_cost_rot_stiefel(Y0_fun, problem_struct_next);
disp(c_Y0_fun)

disp("second RTR output (R_next) -> cost")
c_Rnext = rsom_cost_rot_stiefel(R_next, problem_struct_next);
disp(c_Rnext)

disp("second RTR output projected onto SO(d)^N (R_out) -> cost")
R_out = R_next(1:d, 1:d, 1:N);
c_sci = rsom_cost_rot_stiefel(R_out, problem_struct);
disp(c_sci)

