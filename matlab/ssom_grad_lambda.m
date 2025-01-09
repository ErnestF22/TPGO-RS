function g_lambda = ssom_grad_lambda(X, problem_data)
    edges = problem_data.edges;
    tijs = problem_data.tijs;
    rho = problem_data.rho;

    lambdas = X.lambda;

    g_lambda = zeros(length(lambdas), 1);
    
    num_edges = size(edges, 1);
    for ee = 1:num_edges
        ii = edges(ee, 1);
        jj = edges(ee, 2);
        lambda_e = lambdas(ee);
        tij_e = tijs(:, ee);
        T_i = X.T(:, ii);
        T_j = X.T(:, jj);
        R_i = X.R(:, :, ii);
        a = T_i - T_j;
        b = R_i * tij_e;
        base_part = 2*(b' * b * lambda_e + a' * b);
        relu_part = 0.0;
        if (ssom_relu_argument(lambda_e)>0)
            relu_part = -1.0;
        end
        g_lambda(ee) = base_part + rho * relu_part;
    end
end
