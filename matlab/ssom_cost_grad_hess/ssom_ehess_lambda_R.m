function eh = ssom_ehess_lambda_R(R, Rdot, T, lambdas, problem_data)


% h_lambda_r = zeros(size(h_lambda_lambda));

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
    % lambda_e = lambdas(ee);
    tij = tijs_vec(:, ee);
    T_i = T(:, ii);
    T_j = T(:, jj);
    R_i_dot = Rdot(:, :, ii);
    % R_i = R(:, :, ii);
    
    e_th_elem_half = (tij' * R_i_dot') * (T_i - T_j) ;

    eh(ee) = 2 * e_th_elem_half;
end

end
