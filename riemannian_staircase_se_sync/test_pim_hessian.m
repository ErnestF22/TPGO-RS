close all;
clear;
clc;

resetRands();
%%%
% stream = RandStream('mt19937ar','Seed',0)

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
    'A',A_next,'B',B_next); %should be A, B instead of A_next, B_next


R_manopt_stiefel = stiefelfactory(nrs,d,N);
% fprintf("R_manopt_stiefel size %g\n", R_manopt_stiefel.dim());

step1.M = R_manopt_stiefel; %M = manifold

step1.cost = @(x) som_cost_rot_stiefel(x, problem_struct);
step1.egrad = @(x) som_egrad_rot_stiefel(x, problem_struct);
step1.rgrad = @(x) som_rgrad_rot_stiefel(x, problem_struct);
step1.ehess = @(x,u) som_ehess_rot_stiefel(x,u, problem_struct);
step1.rhess = @(x,u) som_rhess_rot_stiefel(x,u, problem_struct);
%Run Manopt with rand init guess

R_initguess = make_rand_stiefel_3d_array(nrs, d, N);
options.maxiter = 1000;
[R, ~, ~, ~] = trustregions(step1, R_initguess, options);


%%% Run P.I.M. for Hessian

% x = make_rand_stiefel_3d_array(num_rows_stiefel, d, N, array_type);
x = cat_zero_rows_3d_array(R);
u_start = stiefel_randTangentNormVector(x);
stiefel_normalize_han = @(x) x./ (norm(x(:)));
u_start = stiefel_normalize(x, u_start);
rhess_fun_han = @(u) som_rhess_rot_stiefel(x,u,problem_struct_next);


[lambda_pim, v_pim] = pim_function(rhess_fun_han, u_start, stiefel_normalize_han, thresh);
hess_step1 = rhess_fun_han(v_pim);
lambda_pim_prev = lambda_pim;


% v_pim = stiefel_normalize_han(v_pim);
v_pim_prev = v_pim;


disp("lambda_pim found after first iteration of P.I.M.")
disp(lambda_pim)

disp('Difference between lambda*v_max and H(v_max) should be in the order of the tolerance:')
eigencheck_hessian(lambda_pim, v_pim, rhess_fun_han);

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
    eigencheck_hessian(lambda_pim_next, v_pim_next, fun_han_next);

    disp('Checking Eigenvalue shift')
    eigencheck_hessian(lambda_pim_next + mu, v_pim_next, rhess_fun_han);
end

%%%
disp('Checking if highest_norm_eigenval is an eigenval for initial function')
disp('Difference between lambda_min*v_max and H(v_max) should be in the order of the tolerance:')
hess_hne = rhess_fun_han(v_pim);
highest_norm_eigenval = lambda_pim_next + mu;
% disp(norm((highest_norm_eigenval)*v_pim(:) - hess_hne(:),'inf'))
if ~eigencheck_hessian(highest_norm_eigenval, v_pim, rhess_fun_han)
    % scale_factor
    fac_1 = remove_quasi_zeros(highest_norm_eigenval*v_pim(:));
    fac_2 = remove_quasi_zeros(hess_hne(:));
    fac2_1 = fac_2 ./ fac_1;
    fac2_1_nums = fac2_1(~isnan(fac2_1));
    scale_factor = mode(fac2_1_nums);

    disp("Not even after scaling eigenval?")
    eigencheck_hessian(scale_factor * highest_norm_eigenval, v_pim, rhess_fun_han);
    highest_norm_eigenval = scale_factor * highest_norm_eigenval;
end


%Preparing linesearch
step2.M = stiefelfactory(nrs_next,d,N);
step2.cost = @(x) som_cost_rot_stiefel(x, problem_struct_next);
step2.egrad = @(x) som_egrad_rot_stiefel(x, problem_struct_next);
step2.rgrad = @(x) som_rgrad_rot_stiefel(x, problem_struct_next);
step2.ehess = @(x,u) som_ehess_rot_stiefel(x,u, problem_struct_next);
step2.rhess = @(x,u) som_rhess_rot_stiefel(x,u, problem_struct_next);

% alphas = -0.2:0.001:0.2;
% plot_vals = zeros(size(alphas));
% for ii = 1:length(alphas)
%     x_retr_ii = step2.M.retr(x, v_pim, alphas(ii));
%     plot_vals(ii) = step2.cost(x_retr_ii);
% end
% 
% plot(alphas, plot_vals);

alphas = linspace(-0.01,0.01,501); %-0.2:0.01:0.2;
plot_vals = zeros(size(alphas));
plot_vals_taylor = zeros(size(alphas));
for ii = 1:length(alphas)
    x_retr_ii = step2.M.retr(x, v_pim_next, alphas(ii));
    disp("Is x_retr_ii on Stiefel? (Taylor)")
    disp(check_is_on_stiefel(x_retr_ii));
    disp([matStack(x), matStack(x_retr_ii)])
    plot_vals(ii) = step2.cost(x_retr_ii);
    % terms containing the gradient should be zero
    plot_vals_taylor(ii) = step2.cost(x)+...
        alphas(ii)^2/2*sum(stiefel_metric(x,v_pim_next,step2.rhess(x,v_pim_next),'euclidean'));
end

plot(alphas, plot_vals,'b')
hold on
plot(alphas,plot_vals_taylor,'k.');
hold off


% alpha = min(lambdas_moved) + lambdas_max;
alpha_linesearch = 10; %TODO: set this correctly
SDPLRval = 10; %TODO: set this correctly 

% [stepsize, Y0] = linesearch_decrease(step2, ...
%     x, -alpha_linesearch.*v_pim, som_cost_rot_stiefel(x,problem_struct_next));

disp("Now performing linesearch...");
Y0 = linesearch_decrease_hessian(step2, ...
    x, -v_pim, som_cost_rot_stiefel(x,problem_struct_next));

disp("max(abs(x - Y0), [], all)");
disp(max(abs(x - Y0), [], "all"));

%manopt
disp("Reiterating Manopt...");
% R_initguess_next = matUnstack(Y0T, problem_struct_next.sz(1));
R_initguess_next = Y0;
[R_next, R_cost, R_info, R_options] = trustregions(step2, R_initguess_next, options);

% [R_out, T_out] = lowRankLocalization_solution_extractProjection(matStack(multitransp(R_stiefel)) * reshape(T_stiefel, [], N));
% [R_out, T_out] = lowRankLocalization_solution_extractProjection(matStack(multitransp(R_next)));
problem_struct_round_solution = struct("d", d, "n", N);
R_out = round_solution_se_sync(matStackH(R_next), problem_struct_round_solution);


%%% check cost progression

disp("first RTR output (R) -> cost")
c_R = trace(matStack(R)'*problem_struct.L*matStack(R)+matStack(R)'*problem_struct.P) + problem_struct.fixed_cost_term; 
disp(c_R)

disp("first RTR output padded (x) -> cost")
c_x = trace(matStack(x)'*problem_struct_next.L*matStack(x)+matStack(x)'*problem_struct_next.P) + problem_struct_next.fixed_cost_term; 
disp(c_x)

disp("second RTR input (Y0) -> cost")
c_Y0 = trace( ...
    matStack(Y0)'*problem_struct_next.L*matStack(Y0) + ...
    matStack(Y0)'*problem_struct_next.P ) + ...
    problem_struct_next.fixed_cost_term;
disp(c_Y0)

disp("second RTR output (R_next) -> cost")
c_Rnext = trace( ...
    matStack(R_next)'*problem_struct_next.L*matStack(R_next) + ...
    matStack(R_next)'*problem_struct_next.P ) + ...
    problem_struct_next.fixed_cost_term;
disp(c_Rnext)

disp("second RTR output projected onto SO(d)^N (R_out) -> cost")
second_cost_input = matStack(multitransp(matUnstack(R_out')));
c_sci = trace( ...
    second_cost_input'*problem_struct.L*second_cost_input + ...
    second_cost_input'*problem_struct.P ) + ...
    problem_struct.fixed_cost_term;
disp(c_sci)





