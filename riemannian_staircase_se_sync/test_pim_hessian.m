close all;
clear;
clc;

d = 3;
nrs = d;
nrs_next = nrs + 1;
N = 5;
sz = [nrs, d, N];
sz_next = [nrs_next,d,N];

thresh = 1e-5;

array_type = 'double';

L_next = readmatrix("../data/L_stiefel_noisy.csv");
P_next = readmatrix("../data/P_stiefel_noisy.csv");
A_next = readmatrix("../data/A_stiefel_noisy.csv");
B_next = readmatrix("../data/B_stiefel_noisy.csv");
fixed_cost_term = readmatrix("../data/fixed_cost_term.csv");

problem_struct=struct('sz',sz);

problem_struct_next=struct('sz',sz_next, ... %!!
    'L',L_next,'P',P_next, 'fixed_cost_term', fixed_cost_term, ...
    'A',A_next,'B',B_next);


L = from_Lnext_to_L(L_next, problem_struct);
P = from_Pnext_to_P(P_next, problem_struct);
problem_struct=struct('sz',sz, ... %!!
    'L',L,'P',P, 'fixed_cost_term', fixed_cost_term, ...
    'A',A_next,'B',B_next);


R_manopt_stiefel = stiefelfactory(nrs,d,N);
% fprintf("R_manopt_stiefel size %g\n", R_manopt_stiefel.dim());

problem_manopt.M = R_manopt_stiefel; %M = manifold

problem_manopt.cost = @(x) som_cost_rot_stiefel(x, problem_struct);
problem_manopt.egrad = @(x) som_egrad_rot_stiefel(x, problem_struct);
problem_manopt.rgrad = @(x) som_rgrad_rot_stiefel(x, problem_struct);
problem_manopt.ehess = @(x,u) som_ehess_rot_stiefel(x,u, problem_struct);
problem_manopt.rhess = @(x,u) som_rhess_rot_stiefel(x,u, problem_struct);
%Run Manopt with rand init guess

R_initguess = make_rand_stiefel_3d_array(nrs, d, N);
options.maxiter = 100;
[R, R_cost, R_info, R_options] = trustregions(problem_manopt, R_initguess, options);


%%% Run P.I.M. for Hessian

% x = make_rand_stiefel_3d_array(num_rows_stiefel, d, N, array_type);
x = cat_zero_rows_3d_array(R);
u_start = stiefel_randTangentNormVector(x);
stiefel_normalize_han = @(x) x./ (norm(x(:)));
u_start = stiefel_normalize(x, u_start);
fun_han = @(u) som_rhess_rot_stiefel(x,u,problem_struct_next);


[lambda_pim, v_pim] = pim_function(fun_han, u_start, stiefel_normalize_han, thresh);
hess_step1 = fun_han(v_pim);
lambda_pim_prev = lambda_pim;


% v_pim = stiefel_normalize_han(v_pim);
v_pim_prev = v_pim;


disp("lambda_pim found after first iteration of P.I.M.")
disp(lambda_pim)

disp('Difference between lambda*v_max and H(v_max) should be in the order of the tolerance:')
disp(norm(lambda_pim*v_pim(:)- hess_step1(:),'inf'))

%%%%step2
num_iterations = 0;
if lambda_pim>0
    num_iterations = num_iterations + 1;
    fprintf("iteration %g lambda_pim %g\n", num_iterations, lambda_pim);
    
    mu = 1.1 * lambda_pim;

    fun_han_next = @(u) som_rhess_rot_stiefel(x,u,problem_struct_next) - ...
        mu .* u;
            
    %run shifted power iteration
    u_start_new = stiefel_randTangentNormVector(x);
    [lambda_pim_next, v_pim_next] = pim_function(fun_han_next, u_start_new, stiefel_normalize_han, thresh);
    
    disp('Difference between lambda*v_max and H(v_max) should be in the order of the tolerance:')
    hess_next = fun_han_next(v_pim_next);
    disp(norm((lambda_pim_next)*v_pim_next(:) - hess_next(:), 'inf'))

    disp('Checking Eigenvalue shift')
    hess_next2 = som_rhess_rot_stiefel(x,v_pim_next,problem_struct_next);
    disp(norm((lambda_pim+lambda_pim_next)*v_pim_next(:)- hess_next2(:),'inf'))
end

disp("Now performing linesearch...");
%linesearch
step2.M = stiefelfactory(nrs_next,d,N);
step2.cost = @(x) som_cost_rot_stiefel(x, problem_struct_next);
step2.egrad = @(x) som_egrad_rot_stiefel(x, problem_struct_next);
step2.rgrad = @(x) som_rgrad_rot_stiefel(x, problem_struct_next);
step2.ehess = @(x,u) som_ehess_rot_stiefel(x,u, problem_struct_next);
step2.rhess = @(x,u) som_rhess_rot_stiefel(x,u, problem_struct_next);

% alpha = min(lambdas_moved) + lambdas_max;
alpha = 1e-3; %TODO: set this correctly
SDPLRval = 10; %TODO: set this correctly 

[stepsize, Y0T] = linesearch_decrease(step2, ...
    matStack(x), matStack(v_pim), som_cost_rot_stiefel(x,problem_struct_next));

%manopt
disp("Reiterating Manopt...");
% R_initguess_next = make_rand_stiefel_3d_array(nrs_next, d, N);
R_initguess_next = matUnstack(Y0T, problem_struct_next.sz(1));
[R_next, R_cost, R_info, R_options] = trustregions(step2, R_initguess_next, options);

% [R_out, T_out] = lowRankLocalization_solution_extractProjection(matStack(multitransp(R_stiefel)) * reshape(T_stiefel, [], N));
% [R_out, T_out] = lowRankLocalization_solution_extractProjection(matStack(multitransp(R_next)));
problem_struct_round_solution = struct("d", d, "n", N);
R_out = round_solution_se_sync(matStackH(R_next), problem_struct_round_solution);

c = trace(matStack(R)'*problem_struct.L*matStack(R)+matStack(R)'*problem_struct.P) + problem_struct.fixed_cost_term
second_cost_input = matStack(multitransp(matUnstack(R_out'))); 
c = trace( ...
    second_cost_input'*problem_struct.L*second_cost_input + ...
    second_cost_input'*problem_struct.P ) + ...
    problem_struct.fixed_cost_term



