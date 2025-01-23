function h = ssom_ehess_R_T(~, T, Tdot, lambdas, problem_data)
nrs = size(T, 1);
d = size(problem_data.tijs, 1);
N = size(T, 2);

Ph = zeros(nrs, d*N);

tijs_scaled = make_tijs_scaled(lambdas, problem_data.tijs);

idx_col_p = reshape(1:d*N, [], N)';

num_edges = size(problem_data.edges,1);
for e = 1:num_edges
    ii = problem_data.edges(e,1);
    jj = problem_data.edges(e,2);
    Tj_dot = Tdot(:, jj);
    Ti_dot = Tdot(:, ii);
    tij = tijs_scaled(:,e);
    P_e = 2 * (Ti_dot * tij' - Tj_dot * tij');
    Ph(:, idx_col_p(ii, :)) = ...
        Ph(:, idx_col_p(ii, :)) + P_e;
end

h=matUnstackH(Ph,d);

end