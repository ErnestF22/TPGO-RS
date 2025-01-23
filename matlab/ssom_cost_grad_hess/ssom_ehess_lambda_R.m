function eh = ssom_ehess_lambda_R(X, Xdot, problem_data)

R = X.R;
Rdot = Xdot.R;
T = X.T;
lambdas = X.lambda;
% h_lambda_t = zeros(size(h_lambda_lambda));

% x = X.lambda;
% lambdas_dot = Xdot.lambda;
edges = problem_data.edges;
tijs_vec = problem_data.tijs;
% rho = problem_data.rho;

% nrs = problem_data.sz(1);
% % d = problem_data.sz(2);
% N = problem_data.sz(3);

eh = zeros(size(lambdas));

num_edges = size(edges, 1);
for ee = 1:num_edges
    ii = edges(ee, 1);
    jj = edges(ee, 2);
    lambda_e = lambdas(ee);
    tij = tijs_vec(:, ee);
    T_i = T(:, ii);
    T_j = T(:, jj);
    R_i_dot = Rdot(:, :, ii);
    R_i = R(:, :, ii);
    a = T_i - T_j;
    bdot = R_i_dot * tij;
    b = R_i * tij;
    e_th_elem_half = lambda_e * (bdot' * b) + lambda_e * (b' * bdot)  + ...
        a' * bdot;

    eh(ee) = 2 * e_th_elem_half;
end

end
