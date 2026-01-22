function h = ssom_ehess_genproc_relu(X, Xdot, problem_data)
R = X.R;
T = X.T;
lambdas = X.lambda;
Rdot = Xdot.R;
Tdot = Xdot.T;
lambdasdot = Xdot.lambda;

% ehrr = ssom_ehess_R_R(R, Rdot, problem_data); % zero!

% h_r_lambda = zeros(size(R));
hrt = ssom_ehess_R_T(R, T, Tdot, lambdas, problem_data);

h_r_lambda = ssom_ehess_R_lambda(R, T, lambdas, lambdasdot, problem_data);

% h_r_lambda = zeros(size(T));
htr = ssom_ehess_T_R(R, Rdot, T, lambdas, problem_data);

htt = ssom_ehess_T_T(R, T, Tdot, lambdas, problem_data);

h_t_lambda = ssom_ehess_T_lambda(R, T, lambdas, lambdasdot, problem_data);

h_lambda_r = ssom_ehess_lambda_R(R, Rdot, T, lambdas, problem_data);

h_lambda_t = ssom_ehess_lambda_T(R, T, Tdot, lambdas, problem_data);

% h_lambda_lambda = zeros(size(lambda));
h_lambda_lambda = ssom_ehess_lambda_lambda_relu(R, T, lambdas, lambdasdot, problem_data);

% ehR = ehrr + hrt + h_r_lambda;
ehR = hrt + h_r_lambda;
% egR = ssom_egrad_R(R, T, lambdas, problem_data);
h.R = ehR;
h.T = htt + htr + h_t_lambda;
h.lambda = h_lambda_lambda + h_lambda_r + h_lambda_t;
end %rhess genproc

