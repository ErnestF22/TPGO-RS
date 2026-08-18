function test_recompose_eigenvector_from_Hmat

load("data/recompose_stiefel_basis.mat", "X_cat")
load("data/recompose_stiefel_basis.mat", "lambda")
load("data/recompose_stiefel_basis.mat", "v_out_Hmat_tg")
load("data/recompose_stiefel_basis.mat", "rhess_fun_han")

% v_out_Hmat_tg is the column vector corresponding to minimum eigenvalue
v_out_Hmat_struct = recompose_eigenvector_from_Hmat(X_cat, v_out_Hmat_tg);


eigenvalue_check_ok = ssom_eigencheck_hessian_genproc(lambda, v_out_Hmat_struct, rhess_fun_han);

disp("eigenvalue_check_ok")
disp(eigenvalue_check_ok)

end %file function
