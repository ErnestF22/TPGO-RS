function h = ssom_ehess_T_lambda(R, T, ~, lambdas_dot, problem_data)
% h_t_lambda = zeros(size(htr));

edges = problem_data.edges;
% rho = problem_data.rho;

h = zeros(size(T));
N = size(R,3);
num_edges = size(edges, 1);
for e = 1:num_edges
    ii = edges(e,1);
    jj = edges(e,2);
    Ri = R(:,:,ii);
    %         Rj = X.R(:,:,jj);
    %
    BIJ = zeros(N,1);
    BIJ(ii) = 1;
    BIJ(jj) = -1;
    %
    tij = problem_data.tijs(:, e);
    w_ij = BIJ * lambdas_dot(e) * tij' * Ri';
    h = h + 2 * w_ij';
end

end