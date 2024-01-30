function Y0 = rsom_pim_hessian(R, problem_struct_next, thresh)
%RSOM_PIM_HESSIAN Summary of this function goes here
%   Detailed explanation goes here

x = cat_zero_rows_3d_array(R);
u_start = stiefel_randTangentNormVector(x);
stiefel_normalize_han = @(x) x./ (norm(x(:)));
u_start = stiefel_normalize(x, u_start);
rhess_fun_han = @(u) rsom_rhess_rot_stiefel(x,u,problem_struct_next);


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

    fun_han_next = @(u) rsom_rhess_rot_stiefel(x,u,problem_struct_next) - ...
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
nrs_next = problem_struct_next.sz(1);
d = problem_struct_next.sz(2);
N = problem_struct_next.sz(3);
step2.M = stiefelfactory(nrs_next,d,N);
step2.cost = @(x) rsom_cost_rot_stiefel(x, problem_struct_next);
step2.egrad = @(x) rsom_egrad_rot_stiefel(x, problem_struct_next);
step2.grad = @(x) rsom_rgrad_rot_stiefel(x, problem_struct_next);
step2.ehess = @(x,u) rsom_ehess_rot_stiefel(x,u, problem_struct_next);
step2.hess = @(x,u) rsom_rhess_rot_stiefel(x,u, problem_struct_next);

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
%     plot_vals(ii) = step2.cost(x_retr_ii);
%     %Note: gradient is zero
%     plot_vals_taylor(ii) = step2.cost(x)+...
%         alphas(ii)^2/2*sum(stiefel_metric(x,v_pim_next,step2.hess(x,v_pim_next),'canonical'));
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
    x, -v_pim, rsom_cost_rot_stiefel(x,problem_struct_next));

disp("max(abs(x - Y0), [], all)");
disp(max(abs(x - Y0), [], "all"));

% %manopt
% disp("Reiterating Manopt...");
% % R_initguess_next = matUnstack(Y0T, problem_struct_next.sz(1));
% R_initguess_next = Y0;
% [R_next, R_cost, R_info, R_options] = trustregions(step2, R_initguess_next, options);
% 
% % [R_out, T_out] = lowRankLocalization_solution_extractProjection(matStack(multitransp(R_stiefel)) * reshape(T_stiefel, [], N));
% % [R_out, T_out] = lowRankLocalization_solution_extractProjection(matStack(multitransp(R_next)));
% problem_struct_round_solution = struct("d", d, "n", N);
% R_out = round_solution_se_sync(matStackH(R_next), problem_struct_round_solution);


end %file function

