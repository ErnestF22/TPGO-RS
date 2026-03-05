function g_lambda = lsom_egrad_lambda_relu(R, T, lambdas, problem_data)
    edges = problem_data.edges;
    tijs = problem_data.tijs;
    rho = problem_data.rho;

    mu = problem_data.mu;
    y = problem_data.y;
    z = problem_data.z;

    %%
    g_lambda = zeros(length(lambdas), 1);

    num_edges = size(edges, 1);
    for ee = 1:num_edges
        ii = edges(ee, 1);
        jj = edges(ee, 2);
        lambda_e = lambdas(ee);
        tij_e = tijs(:, ee);
        T_i = T(:, ii);
        T_j = T(:, jj);
        R_i = R(:, :, ii);
        base_part = 2*(tij_e'*tij_e * lambda_e + ...
            tij_e' * R_i' * T_i - tij_e' * R_i' * T_j);
        
        if ssom_relu_argument(lambda_e) > 0
            compensation_part = 2 * (lambda_e - 1);
        else
            compensation_part = 0.0;
        end
        
        g_lambda(ee) = base_part + rho * compensation_part;
    end

    lagrange_compensation = -y + mu * (lambdas - z);
    g_lambda = g_lambda + lagrange_compensation;

end
