function g = ssom_egrad_R(R, T, lambdas, problem_data)

g = zeros(size(R));

num_edges = size(problem_data.edges,1);
for e = 1:num_edges
    ii = problem_data.edges(e,1);
    jj = problem_data.edges(e,2);
    Tj = T(:, jj);
    Ti = T(:, ii);
    lambdaij = lambdas(e, :);
    tij = problem_data.tijs(:,e);
    % R_i = R(:,:,ii);
    P_e = 2 * (Ti * lambdaij * tij' - Tj * lambdaij * tij');
    g(:, :, ii) = ...
        g(:, :, ii) + P_e;
end


end %file function