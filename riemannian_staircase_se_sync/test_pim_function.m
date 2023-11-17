d = 3;
num_rows_stiefel = 4;
N = 5;
sz = [num_rows_stiefel, d, N];

array_type = 'double';

L = readmatrix("../data/L_stiefel_noisy.csv");
P = readmatrix("../data/P_stiefel_noisy.csv");
A = readmatrix("../data/A_stiefel_noisy.csv");
B = readmatrix("../data/B_stiefel_noisy.csv");

problem=struct('sz',sz,'L',L,'P',P,'A',A,'B',B);

x = rand_stiefel_tmp(num_rows_stiefel, d, N, array_type);
v_start = stiefel_randTangentNormVector(x);

% [eigvecs, eigvals] = eig(A);

thresh = 1e-5;
fun_han = @(u) som_rhess_rot_stiefel(x,u,problem);
stief_normalize_han = @(x) stiefel_normalize(x);

[lambda_max, v_max] = pim_function(fun_han, v_start, stief_normalize_han, thresh);

disp("v_max")
disp(v_max)

lambda_max = rayleigh_quotient(v_max, A); %Rayleigh quotient

disp("lambda_max")
disp(lambda_max)


% A2 = A - (lambda_max) * eye(size(A));
% 
% [lambda_max2, v_max2] = pim(A2, thresh);
% 
% disp("v_max2")
% disp(v_max2)
% 
% lambda_max2 = rayleigh_quotient(v_max2, A2); %Rayleigh quotient
% 
% disp("lambda_max2")
% disp(lambda_max2)
% 
% disp("lambda_max2 + lambda_max")
% disp(lambda_max2 + lambda_max)


%%%%%%%%%%%%%%%


function rst = rand_stiefel_tmp(n, p, k, array_type)
    rst =  qr_unique(randn(n, p, k, array_type));
end    
    
function U = rand_stiefel_vec_tmp(X, n, p, k, array_type)
    U = stiefel_tangentProj(X, randn(n, p, k, array_type));
    U = U / stiefel_normalize(U(:));
end

