function [Y0, lambda_pim_out, v_pim_out, eigenvalue_check_ok] = ssom_pim_hessian_genproc( ...
    X, problem_struct_next, thresh, num_max_iter)
%SSOM_PIM_HESSIAN_GENPROC Return a new starting point Y0 with lower cost that R
% This is based on a linesearch towards an eigenvector v_pim_out 
% corresponding to negative eigenvalue lambda_pim_out.
% If the Hessian does not have any negative eigenvalue (i.e., the two PIM
% iterations do not find it), simply return Y0 = R and the maximum
% eigenvalue with an associated eigenvector.

if ~exist('thresh', 'var')
    thresh = 1e-6;
end

if ~exist('num_max_iter', 'var')
    num_max_iter = 2000;
end

eigenvalue_check_ok = false;

Rnext = cat_zero_rows_3d_array(X.R);
Tnext = cat_zero_row(X.T);
Xnext.R = Rnext;
Xnext.T = Tnext;
Xnext.lambda = X.lambda;
rhess_fun_han = @(u) ssom_rhess_genproc(Xnext,u,problem_struct_next);

stiefel_normalize_han = @(x) x./ (norm(x(:))); %Note: this is basically eucl_normalize_han

u_start.R = stiefel_randTangentNormVector(Rnext);
u_start.R = stiefel_normalize(Rnext, u_start.R);
u_start.T = rand(size(Tnext));
u_start.T = stiefel_normalize_han(u_start.T);
u_start.lambda = rand(size(Xnext.lambda));
u_start.lambda = stiefel_normalize_han(u_start.lambda);
[lambda_pim, v_pim] = ssom_pim_function_genproc(rhess_fun_han, u_start, stiefel_normalize_han, thresh, num_max_iter);
disp('Difference between lambda*v_max and H(v_max) should be in the order of the tolerance:')
eigenvalue_check_ok = ssom_eigencheck_hessian_genproc(lambda_pim, v_pim, rhess_fun_han);



if lambda_pim>0
%     fprintf("lambda_pim R %g\n", lambda_pim.R);
%     fprintf("lambda_pim T %g\n", lambda_pim.T);
    fprintf("lambda_pim %g\n", lambda_pim);
    
%     lambda_pim = max(lambda_pim.R, lambda_pim.T);

    mu = 1.1 * lambda_pim;

    rhess_shifted_fun_han = ...
        @(u) ssom_rhess_genproc_shifted(Xnext,u,mu,problem_struct_next);
            
    %run shifted power iteration
    u_start_second_iter.R = stiefel_randTangentNormVector(Rnext);
    u_start_second_iter.R = stiefel_normalize(Rnext, u_start_second_iter.R);
    u_start_second_iter.T = rand(size(Tnext));
    u_start_second_iter.T = stiefel_normalize_han(u_start.T);
    u_start_second_iter.lambda = rand(size(Xnext.lambda));
    u_start_second_iter.lambda = stiefel_normalize_han(u_start.lambda);
    [lambda_pim_after_shift, v_pim_after_shift] = ssom_pim_function_genproc( ...
        rhess_shifted_fun_han, u_start_second_iter, stiefel_normalize_han, thresh, num_max_iter);
    
    disp(['Difference between lambda_pim_after_shift*v_pim_after_shift ' ...
        'and H_SH(v_pim_after_shift) should be in the order of the tolerance:'])
    eigenvalue_check_ok = ssom_eigencheck_hessian_genproc(lambda_pim_after_shift, v_pim_after_shift, ...
        rhess_shifted_fun_han);

    disp('Checking Eigenvalue shift:')
    disp(['difference between (lambda_pim_after_shift+mu)*v_pim_after_shift ' ...
        'and H(v_pim_after_shift) should be in the order of the tolerance:'])
    highest_norm_eigenval = lambda_pim_after_shift + mu;
    eigenvalue_check_ok = ssom_eigencheck_hessian_genproc(highest_norm_eigenval, v_pim_after_shift, ...
        rhess_fun_han);
    highest_norm_eigenval = lambda_pim_after_shift + mu;
else
    v_pim_after_shift = v_pim; %variable name in this case is misleading since shift does not happen at all
    highest_norm_eigenval = lambda_pim;
end


%%%
% disp(['Checking if highest_norm_eigenval = lambda_pim_after_shift + mu' ...
%     ' is an eigenval for initial function:'])
% disp(['difference between highest_norm_eigenval*v_pim and H(v_pim) ' ...
%     'should be in the order of the tolerance:'])
 %in case if (lambda_pim>0) FALSE

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
tuple_next.R = euclideanfactory([nrs_next, d, N]);
tuple_next.T = euclideanfactory(nrs_next, N);
num_edges = size(problem_struct_next.edges, 1);
tuple_next.lambda = euclideanfactory(num_edges, 1);
M = productmanifold(tuple_next);
step2.M = M;
step2.sz = [nrs_next, d, N];
step2.cost = @(x) ssom_cost(x, problem_struct_next);
step2.grad = @(x) ssom_egrad(x, problem_struct_next);
step2.hess = @(x, u) ssom_ehess_genproc(x, u, problem_struct_next);

[xRt,dxRt,~,~,ddxRt] = real_geodFun(Xnext.R, v_pim_after_shift.R);
[xTt,dxTt,~,~,ddxTt] = real_geodFun(Xnext.T, v_pim_after_shift.T);
[xLambdat,dxLambdat,~,~,ddxLambdat] = real_geodFun(Xnext.lambda, v_pim_after_shift.lambda);

alphas = linspace(-0.1,0.1,1001); %-0.2:0.01:0.2;
plot_vals = zeros(size(alphas));
plot_vals_taylor = zeros(size(alphas));
for ii = 1:length(alphas)
    % x_retr_ii = step2.M.retr(Xnext, v_pim_after_shift, alphas(ii)); %Manopt retr
%     disp("Is x_retr_ii on Stiefel? (Taylor)")
%     disp(check_is_on_stiefel(x_retr_ii));
%     disp([matStack(x), matStack(x_retr_ii)])

    %formula (5.26) Boumal Intro book;
    % if X is a minimum point, only Hess term (the one without the grads)
    % is non-zero

    x_retr_ii.R = xRt(alphas(ii));
    x_retr_ii.T = xTt(alphas(ii));
    x_retr_ii.lambda = xLambdat(alphas(ii));

    plot_vals(ii) = step2.cost(x_retr_ii);

    %grad
    ssom_rg = ssom_rgrad(Xnext, problem_struct_next);
    rgt_R = alphas(ii)* ...
        sum(stiefel_metric( ...
        Rnext,v_pim_after_shift.R, ...
            ssom_rg.R ,'euclidean'));
    egt_T = alphas(ii)* ...
        sum(stiefel_metric( ...
        Tnext,v_pim_after_shift.T, ...
            ssom_rg.T,'euclidean'));
    egt_lambdas = alphas(ii)* ...
        sum(stiefel_metric( ...
        Xnext.lambda,v_pim_after_shift.lambda, ...
            ssom_rg.lambda,'euclidean'));


    ssom_rhess_var = ssom_rhess_genproc(Xnext, v_pim_after_shift, problem_struct_next);
    pvt_R = alphas(ii)^2/2* ...
        sum(stiefel_metric( ...
        Rnext,v_pim_after_shift.R, ...
            ssom_rhess_var.R ,'euclidean'));
    pvt_T = alphas(ii)^2/2* ...
        sum(stiefel_metric( ...
        Tnext,v_pim_after_shift.T, ...
            ssom_rhess_var.T,'euclidean'));
    pvt_lambdas = alphas(ii)^2/2* ...
        sum(stiefel_metric( ...
        Xnext.lambda,v_pim_after_shift.lambda, ...
            ssom_rhess_var.lambda,'euclidean'));
    
    % OBS. Also following terms are 0 if the curve is a geodesic
    curve_var_R = alphas(ii)^2/2* ...
        sum(stiefel_metric( ...
        Rnext,ssom_rg.R, ...
            ddxRt(0) ,'euclidean'));
    curve_var_T = alphas(ii)^2/2* ...
        sum(stiefel_metric( ...
        Tnext,ssom_rg.T, ...
            ddxTt(0),'euclidean'));
    curve_var_lambdas = alphas(ii)^2/2* ...
        sum(stiefel_metric( ...
        Xnext.lambda,ssom_rg.lambda, ...
            ddxLambdat(0),'euclidean'));
    

    
    plot_vals_taylor(ii) = step2.cost(Xnext) + ...
        rgt_R + egt_T + egt_lambdas + ...
        pvt_R + pvt_T + pvt_lambdas + ...
        curve_var_R + curve_var_T + curve_var_lambdas;
end

figure(100)
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
    Xnext, v_pim_after_shift, ssom_cost(Xnext,problem_struct_next));

% cost_before_ls = ssom_cost(Xnext,problem_struct_next);
% Rnext = Xnext.R(:);
% Tnext = Xnext.T(:);
% lambdasnext = Xnext.lambda(:);
% vpasR = v_pim_after_shift.R(:);
% vpasT = v_pim_after_shift.T(:);
% vpasLambdas = v_pim_after_shift.lambda(:);
% Y0R = Y0.R(:);
% Y0T = Y0.T(:);
% Y0lambda = Y0.lambda(:);
% cost_after_ls = ssom_cost(Y0,problem_struct_next);
% 
% Xnext_vec = [Rnext(:); Tnext(:)];
% Y0_vec = [Y0R(:); Y0T(:)];
% vpas_vec = [vpasR(:); vpasT(:)];
% som_matlab_path = string(getenv("HOME")) + "/workspace/matlab_ws/som/matlab";
% writematrix(cost_before_ls, som_matlab_path + "/data/lsdummy_debug/matlab_cost_before_ls.csv")
% writematrix(Rnext, som_matlab_path + "/data/lsdummy_debug/matlab_Rnext.csv")
% writematrix(Tnext, som_matlab_path + "/data/lsdummy_debug/matlab_Tnext.csv")
% writematrix(vpasR, som_matlab_path + "/data/lsdummy_debug/matlab_vpasR.csv")
% writematrix(vpasT, som_matlab_path + "/data/lsdummy_debug/matlab_vpasT.csv")
% writematrix(Y0R, som_matlab_path + "/data/lsdummy_debug/matlab_Y0R.csv")
% writematrix(Y0T, som_matlab_path + "/data/lsdummy_debug/matlab_Y0T.csv")
% writematrix(cost_after_ls, som_matlab_path + "/data/lsdummy_debug/matlab_cost_after_ls.csv")
% 
% writematrix(Xnext_vec, som_matlab_path + "/data/lsdummy_debug/matlab_Xnext_vec.csv")
% writematrix(vpas_vec, som_matlab_path + "/data/lsdummy_debug/matlab_vpas_vec.csv")
% writematrix(Y0_vec, som_matlab_path + "/data/lsdummy_debug/matlab_Y0_vec.csv")

lambda_pim_out = highest_norm_eigenval;
v_pim_out = v_pim_after_shift;


end %file function

