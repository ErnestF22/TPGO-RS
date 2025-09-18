function test_getVfromVtg_eigencouple
%% not useful for getting eigenvector
load("data/getVfromVtg.mat")


rhess_fun_han = @(u) ssom_rhess_genproc(X_cat,u,problem_data_next);
v_out_Hmat = getVfromVtg(X_cat, v_out_Hmat_tg, lambda_index, problem_data_next);
disp('Difference between lambda*v_max and H(v_max) should be in the order of the tolerance:')
eigenvalue_check_ok = ssom_eigencheck_hessian_genproc(lambda, v_out_Hmat, rhess_fun_han);

disp("eigenvalue_check_ok")
disp(eigenvalue_check_ok)

end %file function





