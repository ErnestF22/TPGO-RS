function h = ssom_ehess_R_lambda(R, T, ~, lambdas_dot, problem_data)

nrs = size(T, 1);
d = size(problem_data.tijs, 1);
N = size(T, 2);

h = zeros(nrs, d, N);

% idx_col_p = reshape(1:d*N, [], N)';

num_edges = size(problem_data.edges,1);
for e = 1:num_edges
    ii = problem_data.edges(e,1);
    jj = problem_data.edges(e,2);
    Tj = T(:, jj);
    Ti = T(:, ii);
    lambdadot_ij = lambdas_dot(e,:);
    tij = problem_data.tijs(:,e);
    R_i = R(:,:,ii);
    P_e = 2 * (Ti * lambdadot_ij * tij' - Tj * lambdadot_ij * tij');
    h(:, :, ii) = ...
        h(:, :, ii) + 2 * P_e - R_i * R_i' * P_e - R_i * P_e' * R_i ;
end

h = 0.5 * h;

end