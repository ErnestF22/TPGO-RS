function h = ssom_ehess_T_R(~, Rdot, ~, lambdas, problem_data)
N = size(Rdot, 3);
nrs = size(Rdot, 1);

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
    lambda_e = lambdas(e, 1);
    % LR = LR + (bij * bij');
    Rdot_i = Rdot(:,:,ii);
    PR_dot_R = PR + lambda_e * (Rdot_i * tij_e * bij');
end

h= 2 * PR_dot_R;
end