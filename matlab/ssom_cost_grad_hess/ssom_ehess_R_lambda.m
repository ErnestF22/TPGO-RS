function h = ssom_ehess_R_lambda(~, T, ~, lambdas_dot, problem_data)

nrs = size(T, 1);
d = size(problem_data.tijs, 1);
N = size(T, 2);

Ph = zeros(nrs, d*N);

tijs_dot_scaled = make_tijs_scaled(lambdas_dot, problem_data.tijs);

idx_col_p = reshape(1:d*N, [], N)';

num_edges = size(problem_data.edges,1);
for e = 1:num_edges
    ii = problem_data.edges(e,1);
    jj = problem_data.edges(e,2);
    T_j = T(:, jj);
    T_i = T(:, ii);
    tij_dot = tijs_dot_scaled(:,e);
    P_e = 2 * (T_i * tij_dot' - T_j * tij_dot');
    Ph(:, idx_col_p(ii, :)) = ...
        Ph(:, idx_col_p(ii, :)) + P_e;
end

h=matUnstackH(Ph,d);

end