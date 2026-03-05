function h = lsom_ehess_lambda_lambda_relu(R, T, lambdas, lambdas_dot, problem_data)
% h_lambda_lambda = zeros(1,1)
% lambdas_dot = Xdot.lambda;
edges = problem_data.edges;
tijs = problem_data.tijs;
% rho = problem_data.rho;

mu = problem_data.mu;

h = zeros(length(lambdas), 1);

num_edges = size(edges, 1);
for ee = 1:num_edges
    % ii = edges(ee, 1);
    % jj = edges(ee, 2);
    % lambda_e = lambdas(ee);
    tij_e = tijs(:, ee);
    % T_i = problem_data.T(:, ii);
    % T_j = problem_data.T(:, jj);
    % R_i = R(:, :, ii);
    % a = T_i - T_j;
    lambda_dot_ee = lambdas_dot(ee,:);


    if ssom_relu_argument(lambdas(ee)) > 0
        compensation_part = 2 * lambda_dot_ee;
    else
        compensation_part = 0.0;
    end

    h(ee) = 2*lambda_dot_ee*(tij_e' * tij_e) + problem_data.rho * compensation_part;
end


h = h + mu * lambdas_dot;

end