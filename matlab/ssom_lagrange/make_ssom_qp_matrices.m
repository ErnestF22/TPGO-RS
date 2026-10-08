function [U, v, A_constr, B_constr] = make_ssom_qp_matrices(R, problem_data)

num_edges = size(problem_data.edges, 1);
edges = problem_data.edges;
N = size(R, 3);

% cost_out = 0.0;

A_full = [];
B_full = [];

for ee = 1:num_edges
    ii = edges(ee, 1);
    jj = edges(ee, 2);

    [A_ee, B_ee] = make_a_b_T_lambda_qp(R, problem_data, ii, jj, ee);

    A_full = [A_full; A_ee];
    B_full = [B_full; B_ee];

end

U = [A_full, B_full];
v = zeros(size(A_full, 1), 1);

A_constr = - [zeros(3*N), zeros(3*N, num_edges); zeros(num_edges, 3*N), eye(num_edges)];
B_constr = [zeros(3*N, 1); -ones(num_edges, 1)];

end