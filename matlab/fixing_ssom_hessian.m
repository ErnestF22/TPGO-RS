function fixing_ssom_hessian

loading = 0;

load('data/test_Hmat_ssom.mat', 'problem_data_next');

if loading
    load('lambda_pim.mat', 'lambda_pim')
    % load('Y_star.mat', 'Y_star')
    % load('v_pim.mat', 'v_pim')
else
    thr = 1e-5;
end

problem_data_next.rho = 0.0;

d = 3;
p = problem_data_next.sz(1);
N = problem_data_next.sz(3);
e = size(problem_data_next.edges, 1);


resetRands(0);
% X.R = eye3d(p,d,N);
X.R = zeros(p,d,N);
resetRands(0);
X.T = rand(p,N);
% X.T = zeros(p,n);
resetRands(0);
X.lambda = rand(e,1);
% X.lambda = zeros(e,1);

X_cat.R = cat_zero_rows_3d_array(X.R);
X_cat.T = cat_zero_row(X.T);
X_cat.lambda = X.lambda;

%%

tuple.R = stiefelfactory(p, d, N);
tuple.T = euclideanfactory(p, N);
tuple.lambda = euclideanfactory(e,1);
M = productmanifold(tuple);
problem_manopt.M = M;
problem_manopt.sz = [p, d, N];
problem_manopt.cost = @(x) ssom_cost(x, problem_data_next);
problem_manopt.grad = @(x) ssom_rgrad(x, problem_data_next);
problem_manopt.hess = @(x, u) ssom_rhess_genproc(x, u, problem_data_next);

v = M.tangent(X, M.rand());
v.R = zeros(size(X.R));

alphas = linspace(-0.01,0.01,501); %-0.2:0.01:0.2;
plot_vals = zeros(size(alphas));
plot_vals_taylor = zeros(size(alphas));
for ii = 1:length(alphas)
    x_retr_ii = problem_manopt.M.retr(X, v, alphas(ii));
%     disp("Is x_retr_ii on Stiefel? (Taylor)")
%     disp(check_is_on_stiefel(x_retr_ii));
%     disp([matStack(x), matStack(x_retr_ii)])
    plot_vals(ii) = problem_manopt.cost(x_retr_ii);
    rhess = ssom_rhess_genproc(X, v, problem_data_next);
    %Note: gradient is zero
    pvt_R = alphas(ii)^2/2* ...
        sum(stiefel_metric( ...
        X.R,v.R, ...
            rhess.R,'canonical'));
    pvt_T = alphas(ii)^2/2* ...
        sum(stiefel_metric( ...
        X.T,v.T, ...
            rhess.T,'euclidean'));
    pvt_lambda = alphas(ii)^2/2* ...
        sum(stiefel_metric( ...
        X.lambda,v.lambda, ...
            rhess.lambda,'euclidean'));
    plot_vals_taylor(ii) = problem_manopt.cost(X) + pvt_R + pvt_T + pvt_lambda;
end

legend
plot(alphas, plot_vals,'b')
hold on
plot(alphas,plot_vals_taylor,'k.');
hold off

%%

if ~loading
    [Y_star, lambda_pim, v_pim] = ssom_pim_hessian_genproc( ...
            X, problem_data_next, thr);
    save('lambda_pim.mat', 'lambda_pim')
    save('Y_star.mat', 'Y_star')
    save('v_pim.mat', 'v_pim')
end
disp("lambda_pim")
disp(lambda_pim)

Hmat_ssom = make_Hmat_ssom(X_cat, problem_data_next);
% disp('Hmat_ssom')
% disp(Hmat_ssom)

[eigvals_Hmat_ssom] = eig(Hmat_ssom);

disp("max(abs(Hmat_ssom - Hmat_ssom'), [], ""all"")")
disp(max(abs(Hmat_ssom - Hmat_ssom'), [], "all"))

% disp("lambda")
% disp(lambda)

disp("min(real(eigvals_Hmat_ssom), [], ""all"")")
disp(min(real(eigvals_Hmat_ssom), [], "all"))



end %file function
