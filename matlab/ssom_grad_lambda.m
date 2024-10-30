function g_lambda = ssom_grad_lambda(x, problem_data)
    edges = problem_data.edges;
    Tijs_vec = problem_data.Tijs;
    rho = problem_data.rho;

    % x = lambdas in this context

    g_lambda = zeros(length(x), 1);
    
    num_edges = size(edges, 1);
    for ee = 1:num_edges
        ii = edges(ee, 1);
        jj = edges(ee, 2);
        lambda_e = x(ee);
        tij_e = Tijs_vec(:, ee);
        T_i = problem_data.T(:, ii);
        T_j = problem_data.T(:, jj);
        R_i = problem_data.R(:, :, ii);
        a = T_i - T_j;
        b = R_i * tij_e;
        base_part = 2*(b' * b * lambda_e + a' * b);
        relu_part = 0.0;
        if (ssom_relu_argument(lambda_e)>0)
            relu_part = ssom_relu_argument(lambda_e);
        end
        g_lambda(ee) = base_part + rho * relu_part;
    end
end