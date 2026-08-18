function cost_out = lsom_cost_relu(X, problem_data)

lambdas = X.lambda;
T = X.T;
R = X.R;

mu = problem_data.mu;
y = problem_data.y;
z = problem_data.z;


edges = problem_data.edges;
tijs = problem_data.tijs;
rho = problem_data.rho;

num_edges = size(edges, 1);

cost_out = 0.0;


% b_log=-a_log/(a_log-1)^2;
for ee = 1:num_edges
    ii = edges(ee, 1);
    jj = edges(ee, 2);
    lambda_e = lambdas(ee);
    tij_e = tijs(:, ee);
    T_i = T(:, ii);
    T_j = T(:, jj);
    R_i = R(:, :, ii);
    a = T_i - T_j;
    b = R_i * tij_e;
    cost_ee = trace(a' * a + 2 * lambda_e * (a' * b) + lambda_e^2 * (b' * b)); 

    scale_compensation_ee = relu_som(ssom_relu_argument(lambda_e));

    cost_out = cost_out + problem_data.rho * scale_compensation_ee + cost_ee;
end

cost_out = cost_out + y'*(vec(z)-vec(lambdas))+0.5 * mu * norm(vec(z)-vec(lambdas))^2;

end %file function