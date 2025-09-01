function h = ssom_ehess_R_T(R, ~, Tdot, lambdas, problem_data)
nrs = size(R, 1);
d = size(problem_data.tijs, 1);
N = size(R, 3);

h = zeros(nrs, d, N);


% idx_col_p = reshape(1:d*N, [], N)';

num_edges = size(problem_data.edges,1);
for e = 1:num_edges
    ii = problem_data.edges(e,1);
    jj = problem_data.edges(e,2);
    Tj_dot = Tdot(:, jj);
    Ti_dot = Tdot(:, ii);
    lambdaij = lambdas(e, :);
    tij = problem_data.tijs(:,e);
    R_i = R(:,:,ii);
    P_e = 2 * (Ti_dot * lambdaij * tij' - Tj_dot * lambdaij * tij');
    h(:, :, ii) = ...
        h(:, :, ii) + 2 * P_e - R_i * R_i' * P_e - R_i * P_e' * R_i ;
end

h = 0.5 * h;

end