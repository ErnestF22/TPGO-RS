function h = lsom_ehess_lambda_lambda(R, T, lambdas, lambdas_dot, problem_data)
% h_lambda_lambda = zeros(1,1)
% lambdas_dot = Xdot.lambda;
edges = problem_data.edges;
tijs = problem_data.tijs;
rho = problem_data.rho;

mu = problem_data.mu;
% y = problem_data.y;
% z = problem_data.z;

a = problem_data.a;
b=-a/(a-1)^2;


h = zeros(length(lambdas), 1);

num_edges = size(edges, 1);
for ee = 1:num_edges
    % ii = edges(ee, 1);
    % jj = edges(ee, 2);
    % lambda_e = lambdas(ee);
    tij_e = tijs(:, ee);
    % T_i = problem_data.T(:, ii);
    % T_j = problem_data.T(:, jj);
    % R_i = R(:, :, ii);
    % a = T_i - T_j;
    lambda_dot_ee = lambdas_dot(ee,:);
    lambda_ee = lambdas(ee,:);

    l = lambda_ee;
    if l <= 1/a
        scale_compensation_ee = 1e+10;
    elseif l<1
        scale_compensation_ee=a/(a*l-1)^2 ...
          +0 ...
          +b;
    else
        scale_compensation_ee=0;
    end

    h(ee) = 2*lambda_dot_ee*(tij_e' * tij_e) + lambda_dot_ee * rho * scale_compensation_ee;
end

h = h + mu * lambdas_dot;

end