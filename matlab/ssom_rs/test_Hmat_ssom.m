function test_Hmat_ssom

close all;

% load('data/test_Hmat_ssom.mat', 'X')
load('data/test_Hmat_ssom.mat', 'problem_data_next')
p = problem_data_next.sz(1);
d = problem_data_next.sz(2);
N = problem_data_next.sz(3);

% p = 2;
% d = 2; 
% N = 2;
% problem_data_next.sz = [p, d, N];
% problem_data_next.edges = [1 2; 2 1];
% num_edges = size(problem_data_next.edges, 1);
% problem_data_next.tijs = 10 * rand(d, num_edges);


num_edges = size(problem_data_next.edges, 1);
problem_data_next.tijs = 10 * rand(d, num_edges);
X.R = make_rand_stiefel_3d_array(p, d, N);
X.T = 50 * rand(p,N);
X.lambda = 5 * rand(num_edges, 1);

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

Hmat_ssom = make_Hmat_ssom_metric(X_cat, problem_data_next);
% disp('Hmat_ssom')
% disp(Hmat_ssom)



% d1 = zeros(size(vectorizeXrtlambdas(X)));
% d1(2,1) = 1;
% d1_struct = convertXtoRTLambdas(d1, p, d, N);
% d2 = zeros(size(vectorizeXrtlambdas(X)));
% d2(1,1) = 1;
% d2_struct = convertXtoRTLambdas(d2, p, d, N);

% lhs = d1' * vectorizeXrtlambdas(ssom_rhess_genproc(X, d2_struct, problem_data_next));
% rhs = d2' * vectorizeXrtlambdas(ssom_rhess_genproc(X, d1_struct, problem_data_next));
% 
% 
% disp("[lhs, Hmat_ssom(2,1)]")
% disp([lhs, Hmat_ssom(2,1)])
% disp("[rhs, Hmat_ssom(1,2)]")
% disp([rhs, Hmat_ssom(1,2)])
% 
% lhs = d1' * vectorizeXrtlambdas(ssom_rhess_genproc(X, d2_struct, problem_data_next));
% rhs = d2' * vectorizeXrtlambdas(ssom_rhess_genproc(X, d1_struct, problem_data_next));

% check_make_H_mat(X_cat, Hmat_ssom, problem_data_next);

[eigvecs_Hmat_ssom, eigvals_Hmat_ssom] = eig(Hmat_ssom);

disp("max(abs(Hmat_ssom - Hmat_ssom'), [], ""all"")")
disp(max(abs(Hmat_ssom - Hmat_ssom'), [], "all"))




[Y0_pim, lambda_pim, v_pim] = ssom_pim_hessian_genproc(X, problem_data_next);


if (is_equal_floats(max(abs(Hmat_ssom - Hmat_ssom'), [], "all"), 0))
    lambda = min(eigvals_Hmat_ssom, [], "all"); % !! HP) Hmat_ssom already symmetric
    
    
    lambda_index = find(lambda == diag(eigvals_Hmat_ssom));
    
    v = eigvecs_Hmat_ssom(:, lambda_index);
    
    disp("v'")
    disp(v')
    
    disp("min(real(eigvals_Hmat_ssom), [], ""all"")")
else
    disp("Hmat_ssom ASYMMETRIC!")
    disp("max(abs(imag(eigvals_Hmat_ssom)), [], ""all"")")
    disp(max(abs(imag(eigvals_Hmat_ssom)), [], "all"))
    lambda = min(real(eigvals_Hmat_ssom), [], "all"); 
end

disp("lambda");
disp(lambda);

disp("lambda_pim")
disp(lambda_pim)

disp("abs(lambda - lambda_pim)")
disp(abs(lambda - lambda_pim))

end %file function



function Hmat = check_make_Hmat(X, Hmat, problem_struct)
% CHECK THE CONSISTENCY IN HMAT BUILDING
Xvec = vectorizeXrtlambdas(X);
vecsz = length(Xvec);
p = size(X.R, 1);
d = size(X.R, 2);
n = size(X.R, 3);
% e = size(problem_struct.edges, 1);
for ii = 1:vecsz
    e_i = zeros(vecsz, 1);
    e_i(ii) = 1;
    % e_i(1:p*d*n) = 1; %(p*d*n + 1:p*d*n + p*n), (p*d*n + p*n + 1:end)
    % e_i(p*d*n + p*n + 1:end) = 1;
    U_e_i = convertXtoRTLambdas(e_i, p, d, n);
    % U_e_i.R = zeros(size(U_e_i.R));
    % U_e_i.lambda = zeros(size(U_e_i.lambda));
    Hgp_e_i = ssom_rhess_genproc(X, U_e_i, problem_struct);
    Hgp_e_i_mat = Hmat * e_i;
    if ~is_equal_floats(vectorizeXrtlambdas(Hgp_e_i), Hgp_e_i_mat)
        disp(is_equal_floats(vectorizeXrtlambdas(Hgp_e_i), Hgp_e_i_mat))
        disp([vectorizeXrtlambdas(Hgp_e_i), Hgp_e_i_mat])
        error("found bug")
    end
end
end

