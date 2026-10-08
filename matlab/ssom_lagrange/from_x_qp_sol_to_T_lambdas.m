function [T,lambda] = from_x_qp_sol_to_T_lambdas(x_opt_qp, N, e)

T = reshape(x_opt_qp(1:3*N), 3, N);
lambda = reshape(x_opt_qp(3*N+1:3*N+e), [], 1);

end
