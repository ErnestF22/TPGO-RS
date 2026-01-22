function h = ssom_ehess_lambda_lambda(R, T, lambdas, lambdas_dot, problem_data)
% h_lambda_lambda = zeros(1,1)
% lambdas_dot = Xdot.lambda;
edges = problem_data.edges;
tijs = problem_data.tijs;
% rho = problem_data.rho;

a = problem_data.a;

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
    lambda_ee = lambdas(ee,:);

    if 1-a*lambda_ee > 0
        compensation_part = (-a*a / (a*a*lambda_ee*lambda_ee - 2 * a * lambda_ee + 1)) + 2 * lambda_dot_ee;
        % compensation_part = a / (1-a*lambda_dot_ee) + 2 * lambda_dot_ee;
    else
        compensation_part = 0.0;
    end

    h(ee) = 2*lambda_dot_ee*(tij_e' * tij_e) + problem_data.rho * compensation_part;
end

end