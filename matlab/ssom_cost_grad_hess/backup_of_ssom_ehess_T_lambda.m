function h = ssom_ehess_T_lambda(R, T, ~, lambdas_dot, problem_data)
% h_t_lambda = zeros(size(htr));

N = size(T, 2);
nrs = size(T, 1);

% LR = zeros(N,N);
PR = zeros(nrs, N);
% BR_const = zeros(d,d);


num_edges = size(problem_data.edges,1);
for e = 1:num_edges
    ii = problem_data.edges(e,1);
    jj = problem_data.edges(e,2);
    bij = zeros(N,1);
    bij(ii, 1) = 1;
    bij(jj, 1) = -1;
    tij_e = problem_data.tijs(:, e);
    lambda_dot_e = lambdas_dot(e, 1);
    % LR = LR + (bij * bij');
    R_i = R(:,:,ii);
    PR_dot_lambda = PR + lambda_dot_e * (R_i * tij_e * bij');
end

h= 2 * PR_dot_lambda;

end