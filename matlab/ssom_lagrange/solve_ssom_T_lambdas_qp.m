function [T, lambdas] = solve_ssom_T_lambdas_qp(R, problem_data)

% --- 1. Define your data matrices U, v ---
% --- 2. Define ONLY Inequality Constraints (A*x <= b) ---

[U, v, A, b] = make_ssom_qp_matrices(R, problem_data);

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
x_optimal = lsqlin(U, v, A, b, Aeq, beq, lb, ub, [], options);

% --- 6. Display Results ---
disp('Optimal solution vector x:');
disp(x_optimal);

% --- 7. Convert the optimal solution to T and lambda values ---
[T, lambdas] = from_x_qp_sol_to_T_lambdas(x_optimal, problem_data.N, problem_data.num_edges);
% disp("T")
% disp(T)
% disp("lambdas")
% disp(lambdas)

end