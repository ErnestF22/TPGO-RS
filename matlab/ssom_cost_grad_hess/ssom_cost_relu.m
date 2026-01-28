function cost_out = ssom_cost_relu(X, problem_data)

lambdas = X.lambda;
T = X.T;
R = X.R;


edges = problem_data.edges;
tijs = problem_data.tijs;
rho = problem_data.rho;

num_edges = size(edges, 1);

cost_out = 0.0;
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
    cost_out = cost_out + cost_ee + rho * scale_compensation_ee * scale_compensation_ee;
    % disp("lambda_e")
    % disp(lambda_e)
    % disp('ee')
    % disp(ee)
    % disp('a')
    % disp(a')
    % disp('b')
    % disp(b')
    % disp('cost_ee')
    % disp(cost_ee)
    % disp('scale_compensation_ee')
    % disp(scale_compensation_ee)
end

end %file function