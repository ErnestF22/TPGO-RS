function g = ssom_egrad_R(~, T, lambdas, problem_data)
nrs = size(T, 1);
d = size(problem_data.tijs, 1);
N = size(T, 2);

P = zeros(nrs, d*N);

tijs_scaled = make_tijs_scaled(lambdas, problem_data.tijs);

idx_col_p = reshape(1:d*N, [], N)';

num_edges = size(problem_data.edges,1);
for e = 1:num_edges
    ii = problem_data.edges(e,1);
    jj = problem_data.edges(e,2); 
    T_j = T(:, jj);
    T_i = T(:, ii);
    tij = tijs_scaled(:,e);
    P_e = 2 * (T_i * tij' - T_j * tij');
    P(:, idx_col_p(ii, :)) = ...
        P(:, idx_col_p(ii, :)) + P_e;
end

g=matUnstackH(P,d);
end