[L_stiefel_prev, P_stiefel_prev] = make_LT_PT_noloops_stiefel(T_globalframe, Tijs_vec, edges, num_rows_stiefel, params);
problem_struct_prev.L = L_stiefel_prev;
problem_struct_prev.P = P_stiefel_prev;
problem_struct_prev.fixed_cost_term = cost_const_term_tij;
problem_struct_prev.num_rows_stiefel = num_rows_stiefel;
problem_struct_prev.sz = [num_rows_stiefel, d, N];
%%%
som_cost_rot_stiefel(R_stiefel, problem_struct_prev)