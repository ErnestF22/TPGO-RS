function h = ssom_rhess_genproc(X, Xdot, problem_data)
R = X.R;
T = X.T;
lambdas = X.lambda;
Rdot = Xdot.R;
Tdot = Xdot.T;
lambdasdot = Xdot.lambda;

hrt = ssom_ehess_R_T(R, T, Tdot, lambdas, problem_data);
htr = ssom_ehess_T_R(R, Rdot, T, lambdas, problem_data);

% h_lambda_lambda = zeros(size(lambda));
h_lambda_lambda = ssom_ehess_lambda_lambda(lambdas, lambdasdot, R, problem_data);

% h_r_lambda = zeros(size(hrt));
h_r_lambda = ssom_ehess_R_lambda(R, T, lambdas, lambdasdot, problem_data);

% h_t_lambda = zeros(size(htr));
h_t_lambda = ssom_ehess_T_lambda(R, T, lambdas, lambdasdot, problem_data);

% h_lambda_r = zeros(size(h_lambda_lambda));
h_lambda_r = ssom_ehess_lambda_R(X, Xdot, problem_data);

% h_lambda_t = zeros(size(h_lambda_lambda));
h_lambda_t = ssom_ehess_lambda_T(R, T, Tdot, lambdas, problem_data);

ehR = ssom_ehess_R_R(R, Rdot, problem_data) + hrt + h_r_lambda;
egR = ssom_egrad_R(R, T, lambdas, problem_data);
h.R = manopt_stiefel_ehess2rhess(R, egR, ehR, Rdot);
h.T = ssom_ehess_T_T(R, T, Tdot, lambdas, problem_data) + htr + h_t_lambda;
h.lambda = h_lambda_lambda + h_lambda_r + h_lambda_t;
end %rhess genproc