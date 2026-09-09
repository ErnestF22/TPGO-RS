function test_solve_ssom_T_lambdas_qp

sigma_noise = 0.0;

N = 5;
mindeg = 2;
testdata = testNetwork_params(3, N, 'banded', mindeg, 0.0, sigma_noise);
testdata.mindeg = mindeg;
testdata.sigma = sigma_noise;
testdata.rho = 0; % compensation part
testdata.a = 2.0; % log scale-compensation param

num_edges = size(testdata.E, 1);


X.R = G2R(testdata.gi);
X.T = G2T(testdata.gi);
X.lambda = testdata.lambdaij;

% resetRands(posixtime(datetime("now")));
% 
% X.R = rand(size(X.R));
% X.T = rand(size(X.T));
% X.lambda = rand(size(X.lambda));


problem_data.tijs = G2T(testdata.gij);
problem_data.tijs = rand(size(problem_data.tijs));
problem_data.edges = testdata.E;
problem_data.rho = 0.0;
problem_data.a = testdata.a;

% --- 1. Define your data matrices U, v ---
% --- 2. Define ONLY Inequality Constraints (A*x <= b) ---

[U, v, A, b] = make_ssom_qp_matrices(X.R, problem_data);

% --- 3. Define Variable Bounds (Optional Inequality) ---
% If your variables must be non-negative (x >= 0), set lb = 0.
% If there are no bounds, leave them as empty arrays [].
lb = []; 
ub = [];

% --- 4. Leave Equality Fields Empty ---
Aeq = [];
beq = [];

% --- 5. Solve the problem ---
options = optimoptions('lsqlin', 'Display', 'final', 'Algorithm', 'interior-point');
x_optimal = lsqlin(U, v, -A, b, Aeq, beq, lb, ub, [], options);

% --- 6. Display Results ---
disp('Optimal solution vector x:');
disp(x_optimal);

% --- 7. Convert the optimal solution to T and lambda values ---
[T, lambda] = from_x_qp_sol_to_T_lambdas(x_optimal, N, num_edges);
disp("T")
disp(T)
disp("lambda")
disp(lambda)

end %file function

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

A_constr = [zeros(3*N), zeros(3*N, num_edges); zeros(num_edges, 3*N), eye(num_edges)];
B_constr = [zeros(3*N, 1); -ones(num_edges, 1)];

end

function [T,lambda] = from_x_qp_sol_to_T_lambdas(x_opt_qp, N, e)

T = reshape(x_opt_qp(1:3*N), 3, N);
lambda = reshape(x_opt_qp(3*N+1:3*N+e), [], 1);

end
