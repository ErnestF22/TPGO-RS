function h = ssom_ehess_lambda_T(R, ~,  Tdot, lambdas, problem_data)

% h_lambda_r = zeros(size(h_lambda_lambda));

% x = X.lambdas;
% lambdas_dot = Xdot.lambdas;
edges = problem_data.edges;
tijs_vec = problem_data.tijs;
% rho = problem_data.rho;

h = zeros(size(lambdas));

num_edges = size(edges, 1);
for ee = 1:num_edges
    ii = edges(ee, 1);
    jj = edges(ee, 2);
    % lambda_e = x(ee);
    tij = tijs_vec(:, ee);
    % T_i = T(:,ii);
    % T_j = T(:,jj);
    T_i_dot = Tdot(:, ii);
    T_j_dot = Tdot(:, jj);
    % a = T_i - T_j;
    R_i = R(:, :, ii);
    b = R_i * tij;
    adot = T_i_dot - T_j_dot;
    e_th_elem = 2 * adot' * b;
    h(ee) = e_th_elem;
end

end