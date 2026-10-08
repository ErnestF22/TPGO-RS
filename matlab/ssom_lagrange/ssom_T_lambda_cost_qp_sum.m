function cost_out = ssom_T_lambda_cost_qp_sum(X, problem_data)
num_edges = size(problem_data.edges, 1);
edges = problem_data.edges;
R = X.R;
T = X.T;
lambdas = X.lambda;
cost_out = 0.0;

for ee = 1:num_edges
    ii = edges(ee, 1);
    jj = edges(ee, 2);

    [A_ee, B_ee] = make_a_b_T_lambda_qp(R, problem_data, ii, jj, ee);

    cost_increment = norm([A_ee, B_ee] * [T(:); lambdas(:)]);
    cost_out = cost_out + cost_increment*cost_increment;

end

end