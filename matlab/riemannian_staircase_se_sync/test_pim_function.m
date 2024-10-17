d = 3;
num_rows_stiefel = 4;
N = 5;
sz = [num_rows_stiefel, d, N];

array_type = 'double';

L = readmatrix("../data/L_stiefel_noisy.csv");
P = readmatrix("../data/P_stiefel_noisy.csv");
A = readmatrix("../data/A_stiefel_noisy.csv");
B = readmatrix("../data/B_stiefel_noisy.csv");

problem_struct=struct('sz',sz,'L',L,'P',P,'A',A,'B',B);

x = make_rand_stiefel_3d_array(num_rows_stiefel, d, N);
v_start = stiefel_randTangentNormVector(x);

% [eigvecs, eigvals] = eig(A);

thresh = 1e-5;
fun_han = @(u) som_rhess_rot_stiefel(x,u,problem_struct);
% stief_normalize_han = @(x) stiefel_normalize(x);
stiefel_normalize_han = @(x) x./ (norm(x(:)));

[lambda_max, v_max] = pim_function(fun_han, v_start, stiefel_normalize_han, thresh);

disp("v_max")
disp(v_max)

disp("lambda_max")
disp(lambda_max)

%check if lambda_pim is actually an eigenvalue:
disp("is lambda_max an eigenvalue?")
disp([lambda_max * v_max, som_rhess_rot_stiefel(x, v_max, problem_struct)])


x_next = cat_zero_rows_3d_array(x);
v_start_next = stiefel_randTangentNormVector(x_next);
% v_max_next = cat_zero_rows_3d_array(v_max);
problem_next_struct = problem_struct;
problem_next_struct.sz(1) = problem_next_struct.sz(1) + 1;
% stief_normalize_han = @(x) stiefel_normalize(x);
L_next = from_L_to_L_next(problem_struct.L, problem_next_struct);
P_next = from_P_to_P_next(problem_struct.P, problem_next_struct);
stiefel_normalize_han = @(x) x./ (norm(x(:)));
problem_next_struct.L = L_next;
problem_next_struct.P = P_next;

fun_han_next = @(u) som_rhess_rot_stiefel(x_next,u,problem_next_struct);
[lambda_pim_next, v_max_next] = pim_function(fun_han_next, v_start_next, stiefel_normalize_han, thresh);



%%%%%%%%%%%%%%%

    
function U = rand_stiefel_vec_tmp(X, n, p, k, array_type)
    U = stiefel_tangentProj(X, randn(n, p, k, array_type));
    U = U / stiefel_normalize(U(:));
end

