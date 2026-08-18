function g_lambda = lsom_egrad_lambda(R, T, lambdas, problem_data)
    edges = problem_data.edges;
    tijs = problem_data.tijs;

    mu = problem_data.mu;
    y = problem_data.y;
    z = problem_data.z;

    g_lambda = zeros(length(lambdas), 1);

    rho = problem_data.rho;
    a = problem_data.a;
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
        
        l = lambda_e;
        if l <= 1/a
            scale_compensation_ee = 1e+10;
        elseif l<1
            scale_compensation_ee=-1/(a*l-1)...
                +1/(a-1)...
                +b*(l-1);
        else
            scale_compensation_ee=0;
        end

        % y_ee = y(ee,1);
        % z_ee = z(ee,1);        

        g_lambda(ee) = base_part + rho * scale_compensation_ee;
    end

    lagrange_compensation = -y + mu * (lambdas - z);
    g_lambda = g_lambda + lagrange_compensation;

end %file function
