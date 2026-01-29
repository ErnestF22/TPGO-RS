function g_lambda = ssom_egrad_lambda(R, T, lambdas, problem_data)
    edges = problem_data.edges;
    tijs = problem_data.tijs;
    rho = problem_data.rho;
    a = problem_data.a;


    g_lambda = zeros(length(lambdas), 1);

    b=-a/(a-1)^2;
    
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
        
        % if 1-a*lambda_e > 0
        %     scale_compensation_ee = a / (1-a*lambda_e) - a + 2*lambda_e;
        % else
        %     scale_compensation_ee = 0.0;
        % end

        l = lambda_e;
        if l<=1
            scale_compensation_ee=-1/(a*l-1)...
                +1/(a-1)...
                +b*(l-1);
        else
            scale_compensation_ee=0;
        end
        
        g_lambda(ee) = base_part + rho * scale_compensation_ee;
    end
end
