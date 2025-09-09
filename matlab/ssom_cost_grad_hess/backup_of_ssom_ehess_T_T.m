function h = ssom_ehess_T_T(R, T, Tdot, lambdas, problem_data)


N = size(T, 2);
% nrs = size(T, 1);

LR = zeros(N,N);
% PR = zeros(N,nrs);
% BR_const = zeros(d,d);


num_edges = size(problem_data.edges,1);
for e = 1:num_edges
    ii = problem_data.edges(e,1);
    jj = problem_data.edges(e,2);
    bij = zeros(N,1);
    bij(ii, 1) = 1;
    bij(jj, 1) = -1;
    % tij = problem_data.tijs(:, e);
    % lambda_e = lambdas(e, 1);
    LR = LR + (bij * bij');
    % Ri = R(:,:,ii);
    % PR = PR + 2 * (bij * lambda_e * tij' * Ri');
end

h = Tdot*(LR' + LR);
end
