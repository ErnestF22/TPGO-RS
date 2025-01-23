function h = ssom_ehess_lambda_lambda(lambdas, lambdasdot, R, problem_data)
% h_lambda_lambda = zeros(1,1)
% lambdas_dot = Xdot.lambda;
edges = problem_data.edges;
tijs = problem_data.tijs;
% rho = problem_data.rho;

h = zeros(length(lambdas), 1);

num_edges = size(edges, 1);
for ee = 1:num_edges
    ii = edges(ee, 1);
    % jj = edges(ee, 2);
    % lambda_e = lambdas(ee);
    tij_e = tijs(:, ee);
    % T_i = problem_data.T(:, ii);
    % T_j = problem_data.T(:, jj);
    R_i = R(:, :, ii);
    % a = T_i - T_j;
    b = R_i * tij_e;
    h(ee) = 2*(b' * b) * lambdasdot(ee);
end

end