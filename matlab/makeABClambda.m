function [aL, bL, cL] = makeABClambda(X, problem_data)
lambdas = X.lambda;
T = X.T;
R = X.R;


edges = problem_data.edges;
tijs_vec = problem_data.tijs;
% rho = problem_data.rho;

num_edges = size(edges, 1);
aL = zeros(size(lambdas));
bL = zeros(size(lambdas));
cL = zeros(size(lambdas));

% cost_out = 0.0;
for ee = 1:num_edges
    ii = edges(ee, 1);
    jj = edges(ee, 2);
    % lambda_e = lambdas(ee);
    tij_e = tijs_vec(:, ee);
    T_i = T(:, ii);
    T_j = T(:, jj);
    R_i = R(:, :, ii);
    a = T_i - T_j;
    b = R_i * tij_e;
    % cost_lambda_0_ee = trace(aL' * aL + 2 * lambda_e * (aL' * b) + lambda_e^2 * (b' * b));
    % cost_relu_ee = relu_som(ssom_relu_argument(lambda_e));
    % cost_out = cost_out + cost_lambda_0_ee + rho * cost_relu_ee;
    aL(ee) = trace(a' * a);
    bL(ee) = 2 * trace(a' * b);
    cL(ee) = trace(b' * b);
end
end %function