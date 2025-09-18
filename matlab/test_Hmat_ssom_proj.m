function [lambda, lambda_pim, v_out_Hmat, v_out_pim, pim_eigenvalue_check_ok, imag_eigenvalues] = test_Hmat_ssom_proj

load('data/test_Hmat_ssom.mat', 'problem_data_next')
problem_data = problem_data_next; % !!
p = problem_data.sz(1);
d = problem_data.sz(2);
N = problem_data.sz(3);
problem_data_next.sz(1) = problem_data.sz(1) + 1;



num_edges = size(problem_data_next.edges, 1);
problem_data_next.tijs = 10 * rand(d, num_edges);
X.R = make_rand_stiefel_3d_array(p, d, N);
X.T = 50 * rand(p,N);
X.lambda = 5 * rand(num_edges, 1);

% stb = stiefel_tangentBasis(X.R(:,:,1));

problem_data_next.rho = 0;


Xvec = vectorizeXrtlambdas(X);
X2 = convertXtoRTLambdas(Xvec, p, d, N);

% checking vectorizeXrtlambdas() and convertXtoRTLambdas()
disp("checking vectorizeXrtlambdas() and convertXtoRTLambdas()")
disp("is_equal_floats(X.R, X2.R)")
disp(is_equal_floats(X.R, X2.R))
disp("is_equal_floats(X.T, X2.T)")
disp(is_equal_floats(X.T, X2.T))
disp("is_equal_floats(X.lambda, X2.lambda)")
disp(is_equal_floats(X.lambda, X2.lambda))

% problem_data_next.rho = -10.0;

% R = X.R;
% T = X.T;
% Lambda = X.lambda;

X_cat.R = cat_zero_rows_3d_array(X.R);
X_cat.T = cat_zero_row(X.T);
X_cat.lambda = X.lambda;

% X_cat.R = zeros(size(cat_zero_rows_3d_array(X.R)));
% X_cat.T = rand(size(cat_zero_row(X.T)));
% X_cat.lambda = rand(size(X.lambda));

Hmat_ssom = make_Hmat_ssom_proj(X_cat, problem_data_next);
% disp('Hmat_ssom')
% disp(Hmat_ssom)


[eigvecs_Hmat_ssom, eigvals_Hmat_ssom] = eig(Hmat_ssom);

disp("max(abs(Hmat_ssom - Hmat_ssom'), [], ""all"")")
disp(max(abs(Hmat_ssom - Hmat_ssom'), [], "all"))




[Y0_pim, lambda_pim, v_out_pim, pim_eigenvalue_check_ok] = ...
    ssom_pim_hessian_genproc(X, problem_data_next, 1e-6, 5000);

imag_eigenvalues = false;

if (is_equal_floats(max(abs(Hmat_ssom - Hmat_ssom'), [], "all"), 0))
    lambda = min(eigvals_Hmat_ssom, [], "all"); % !! HP) Hmat_ssom already symmetric
    
    
    lambda_index = find(lambda == diag(eigvals_Hmat_ssom));
    
    v_out_Hmat = eigvecs_Hmat_ssom(:, lambda_index);
    
    % disp("v'")
    % disp(v')
    
    disp("min(real(eigvals_Hmat_ssom), [], ""all"")")
else
    disp("Hmat_ssom ASYMMETRIC!")
    disp("max(abs(imag(eigvals_Hmat_ssom)), [], ""all"")")
    disp(max(abs(imag(eigvals_Hmat_ssom)), [], "all"))
    if max(abs(imag(eigvals_Hmat_ssom)), [], "all") > 1e-5
        imag_eigenvalues = true;
    end
    lambda = min(real(eigvals_Hmat_ssom), [], "all"); 
    lambda_index = find(lambda == diag(real(eigvals_Hmat_ssom)));
    
    v_out_Hmat_tg = real(eigvecs_Hmat_ssom(:, lambda_index));
end

disp("lambda");
disp(lambda);

disp("lambda_pim")
disp(lambda_pim)

%%
rhess_fun_han = @(u) ssom_rhess_genproc(X_cat,u,problem_data_next);
v_out_Hmat = getVfromVtg(X_cat, v_out_Hmat_tg, lambda_index, problem_data_next);
disp('Difference between lambda*v_max and H(v_max) should be in the order of the tolerance:')
eigenvalue_check_ok = ssom_eigencheck_hessian_genproc(lambda, v_out_Hmat, rhess_fun_han);


end %file function
