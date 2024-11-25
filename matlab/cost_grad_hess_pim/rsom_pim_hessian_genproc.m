function [Y0, lambda_pim_out, v_pim_out] = rsom_pim_hessian_genproc( ...
    X, problem_struct_next, thresh)
%RSOM_PIM_HESSIAN Return a new starting point Y0 with lower cost that R
% This is based on a linesearch towards an eigenvector v_pim_out 
% corresponding to negative eigenvalue lambda_pim_out.
% If the Hessian does not have any negative eigenvalue (i.e., the two PIM
% iterations do not find it), simply return Y0 = R and the maximum
% eigenvalue with an associated eigenvector.

if ~exist('thresh', 'var')
    thresh = 1e-6;
end

Rnext = cat_zero_rows_3d_array(X.R);
Tnext = cat_zero_row(X.T);
Xnext.R = Rnext;
Xnext.T = Tnext;
rhess_fun_han = @(u) hess_genproc(Xnext,u,problem_struct_next);

stiefel_normalize_han = @(x) x./ (norm(x(:))); %Note: this is basically eucl_normalize_han

u_start.R = stiefel_randTangentNormVector(Rnext);
u_start.R = stiefel_normalize(Rnext, u_start.R);
u_start.T = rand(size(Tnext));
u_start.T = stiefel_normalize_han(u_start.T);
[lambda_pim, v_pim] = pim_function_genproc(rhess_fun_han, u_start, stiefel_normalize_han, thresh);
disp('Difference between lambda*v_max and H(v_max) should be in the order of the tolerance:')
eigencheck_hessian_genproc(lambda_pim, v_pim, rhess_fun_han);



if lambda_pim>0
%     fprintf("lambda_pim R %g\n", lambda_pim.R);
%     fprintf("lambda_pim T %g\n", lambda_pim.T);
    fprintf("lambda_pim %g\n", lambda_pim);
    
%     lambda_pim = max(lambda_pim.R, lambda_pim.T);

    mu = 1.1 * lambda_pim;

    rhess_shifted_fun_han = ...
        @(u) hess_genproc_shifted(Xnext,u,mu,problem_struct_next);
            
    %run shifted power iteration
    u_start_second_iter.R = stiefel_randTangentNormVector(Rnext);
    u_start_second_iter.R = stiefel_normalize(Rnext, u_start_second_iter.R);
    u_start_second_iter.T = rand(size(Tnext));
    u_start_second_iter.T = stiefel_normalize_han(u_start.T);
    [lambda_pim_after_shift, v_pim_after_shift] = pim_function_genproc( ...
        rhess_shifted_fun_han, u_start_second_iter, stiefel_normalize_han, thresh);
    
    disp(['Difference between lambda_pim_after_shift*v_pim_after_shift ' ...
        'and H_SH(v_pim_after_shift) should be in the order of the tolerance:'])
    eigencheck_hessian_genproc(lambda_pim_after_shift, v_pim_after_shift, ...
        rhess_shifted_fun_han);

    disp('Checking Eigenvalue shift:')
    disp(['difference between (lambda_pim_after_shift+mu)*v_pim_after_shift ' ...
        'and H(v_pim_after_shift) should be in the order of the tolerance:'])
    highest_norm_eigenval = lambda_pim_after_shift + mu;
    eigencheck_hessian_genproc(highest_norm_eigenval, v_pim_after_shift, ...
        rhess_fun_han);
end
%TODO! if lambda_pim already < 0

%%%
% disp(['Checking if highest_norm_eigenval = lambda_pim_after_shift + mu' ...
%     ' is an eigenval for initial function:'])
% disp(['difference between highest_norm_eigenval*v_pim and H(v_pim) ' ...
%     'should be in the order of the tolerance:'])
highest_norm_eigenval = lambda_pim_after_shift + mu; %in case if (lambda_pim>0) FALSE

% eigencheck_hessian_genproc(highest_norm_eigenval, v_pim, rhess_fun_han);
%%% scaling eigenvalue
% if ~eigencheck_hessian(highest_norm_eigenval, v_pim, rhess_fun_han)
%     % scale_factor
%     fac_1 = remove_quasi_zeros(highest_norm_eigenval*v_pim(:));
%     hess_hne = rhess_fun_han(v_pim);
%     fac_2 = remove_quasi_zeros(hess_hne(:));
%     fac2_1 = fac_2 ./ fac_1;
%     fac2_1_nums = fac2_1(~isnan(fac2_1));
%     fac2_1_finite = fac2_1_nums(isfinite(fac2_1_nums));
%     scale_factor = mode(fac2_1_finite);
% 
%     disp("Not even after scaling eigenval?")
%     eigencheck_hessian(scale_factor * highest_norm_eigenval, v_pim, rhess_fun_han);
%     highest_norm_eigenval = scale_factor * highest_norm_eigenval;
% end


%Preparing linesearch
nrs_next = problem_struct_next.sz(1);
d = problem_struct_next.sz(2);
N = problem_struct_next.sz(3);
tuple_next.R = stiefelfactory(nrs_next, d, N);
tuple_next.T = euclideanfactory(nrs_next, N);
M = productmanifold(tuple_next);
step2.M = M;
step2.sz = [nrs_next, d, N];
step2.cost = @(x) cost_genproc(x, problem_struct_next);
step2.grad = @(x) grad_genproc(x, problem_struct_next);
step2.hess = @(x, u) hess_genproc(x, u, problem_struct_next);


alphas = linspace(-0.01,0.01,501); %-0.2:0.01:0.2;
plot_vals = zeros(size(alphas));
plot_vals_taylor = zeros(size(alphas));
for ii = 1:length(alphas)
    x_retr_ii = step2.M.retr(Xnext, v_pim_after_shift, alphas(ii));
%     disp("Is x_retr_ii on Stiefel? (Taylor)")
%     disp(check_is_on_stiefel(x_retr_ii));
%     disp([matStack(x), matStack(x_retr_ii)])
    plot_vals(ii) = step2.cost(x_retr_ii);
    %Note: gradient is zero
    pvt_R = alphas(ii)^2/2* ...
        sum(stiefel_metric( ...
        Rnext,v_pim_after_shift.R, ...
            hess_genproc_R(Xnext, v_pim_after_shift, problem_struct_next),'canonical'));
    pvt_T = alphas(ii)^2/2* ...
        sum(stiefel_metric( ...
        Tnext,v_pim_after_shift.T, ...
            hess_genproc_T(Xnext, v_pim_after_shift, problem_struct_next),'euclidean'));
    plot_vals_taylor(ii) = step2.cost(Xnext) + pvt_R + pvt_T;
end

plot(alphas, plot_vals,'b')
hold on
plot(alphas,plot_vals_taylor,'k.');
hold off


% alpha = min(lambdas_moved) + lambdas_max;
% alpha_linesearch = 10; %TODO: set this correctly
% SDPLRval = 10; %TODO: set this correctly 

disp("Now performing linesearch...");
%Note: first output param of linesearch() would be "stepsize"
[~, Y0] = linesearch_decrease(step2, ...
    Xnext, v_pim_after_shift, cost_genproc(Xnext,problem_struct_next));

lambda_pim_out = highest_norm_eigenval;
v_pim_out = v_pim_after_shift;


end %file function

