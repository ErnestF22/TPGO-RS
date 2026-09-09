function cost_out = ssom_T_lambda_cost_qp_matricial(X, problem_data)
num_edges = size(problem_data.edges, 1);
edges = problem_data.edges;
R = X.R;
T = X.T;
lambdas = X.lambda;
% cost_out = 0.0;

A_ee_full = [];
B_ee_full = [];

for ee = 1:num_edges
    ii = edges(ee, 1);
    jj = edges(ee, 2);

    [A_ee, B_ee] = make_a_b_T_lambda_qp(R, problem_data, ii, jj, ee);

    A_ee_full = [A_ee_full; A_ee];
    B_ee_full = [B_ee_full; B_ee];

end


cost_out = norm([A_ee_full, B_ee_full] * [T(:); lambdas(:)]); % Accumulate the cost

cost_out = cost_out * cost_out;

end