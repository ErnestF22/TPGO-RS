function h = ssom_ehess_genproc(X, Xdot, problem_data)
R = X.R;
T = X.T;
lambdas = X.lambda;
Rdot = Xdot.R;
Tdot = Xdot.T;
lambdasdot = Xdot.lambda;

% hrt = ssom_ehess_R_T(R, T, Tdot, lambdas, problem_data);
hrt = computeHrt(R, T, Tdot, lambdas, problem_data);
htr = ssom_ehess_T_R(R, Rdot, T, lambdas, problem_data);

% h_lambda_lambda = zeros(size(lambda));
h_lambda_lambda = ssom_ehess_lambda_lambda(R, T, lambdas, lambdasdot, problem_data);

% h_r_lambda = zeros(size(hrt));
% h_r_lambda = ssom_ehess_R_lambda(R, T, lambdas, lambdasdot, problem_data);
h_r_lambda = computeHrlambda(R, T, lambdas, lambdasdot, problem_data);

% h_t_lambda = zeros(size(htr));
h_t_lambda = ssom_ehess_T_lambda(R, T, lambdas, lambdasdot, problem_data);

% h_lambda_r = zeros(size(h_lambda_lambda));
h_lambda_r = ssom_ehess_lambda_R(R, Rdot, T, lambdas, problem_data);

% h_lambda_t = zeros(size(h_lambda_lambda));
h_lambda_t = ssom_ehess_lambda_T(R, T, Tdot, lambdas, problem_data);

ehR = ssom_ehess_R_R(R, Rdot, problem_data) + hrt + h_r_lambda;
% egR = ssom_egrad_R(R, T, lambdas, problem_data);
h.R = ehR;
h.T = ssom_ehess_T_T(R, T, Tdot, lambdas, problem_data) + htr + h_t_lambda;
h.lambda = h_lambda_lambda + h_lambda_r + h_lambda_t;
end %rhess genproc

function hrt = computeHrt(R, T, Tdot, lambdas, problem_data)
    nrs = size(R, 1);
    d = size(problem_data.tijs, 1);
    N = size(R, 3);
    
    hrt = zeros(nrs, d, N);
        
    num_edges = size(problem_data.edges,1);
    for e = 1:num_edges
        ii = problem_data.edges(e,1);
        jj = problem_data.edges(e,2);
        Tj_dot = Tdot(:, jj);
        Ti_dot = Tdot(:, ii);
        lambdaij = lambdas(e, :);
        tij = problem_data.tijs(:,e);
        % R_i = R(:,:,ii);
        P_e = 2 * (Ti_dot * lambdaij * tij' - Tj_dot * lambdaij * tij');
        hrt(:, :, ii) = ...
            hrt(:, :, ii) + P_e;
    end
    
end

function hrlambda = computeHrlambda(R, T, lambdas, lambdasdot, problem_data)

    nrs = size(R, 1);
    d = size(problem_data.tijs, 1);
    N = size(R, 3);
    
    hrlambda = zeros(nrs, d, N);    
        
    num_edges = size(problem_data.edges,1);
    for e = 1:num_edges
        ii = problem_data.edges(e,1);
        jj = problem_data.edges(e,2);
        Tj = T(:, jj);
        Ti = T(:, ii);
        lambdaij_dot = lambdasdot(e, :);
        tij = problem_data.tijs(:,e);
        % R_i = R(:,:,ii);
        P_e = 2 * (Ti * lambdaij_dot * tij' - Tj * lambdaij_dot * tij');
        hrlambda(:, :, ii) = ...
            hrlambda(:, :, ii) + P_e;
    end
    
end