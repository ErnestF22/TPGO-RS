function cost_out = ssom_cost_norms(X, problem_data)

lambdas = X.lambda;
T = X.T;
R = X.R;


edges = problem_data.edges;
Tijs_vec = problem_data.Tijs;
rho = problem_data.rho;

num_edges = size(edges, 1);

cost_out = 0.0;
for ee = 1:num_edges
    ii = edges(ee,1);
    jj = edges(ee,2); 
    R_i = R(:,:,ii);
    T_j = T(:, jj);
    T_i = T(:, ii);
    lambda_e = lambdas(ee);
    tij_e = Tijs_vec(:, ee);
    cost_e = norm(R_i * lambda_e * tij_e - T_j + T_i);
    cost_relu_ee = relu_som(ssom_relu_argument(lambda_e));
    cost_out = cost_out + cost_e^2 + rho * cost_relu_ee; % cost_e squared!
end

end %file function
