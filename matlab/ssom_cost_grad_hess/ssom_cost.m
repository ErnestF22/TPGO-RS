function cost_out = ssom_cost(X, problem_data)

lambdas = X.lambda;
T = X.T;
R = X.R;


edges = problem_data.edges;
tijs = problem_data.tijs;
rho = problem_data.rho;
a_log = problem_data.a;

num_edges = size(edges, 1);

cost_out = 0.0;

b_log=-a_log/(a_log-1)^2;
for ee = 1:num_edges
    ii = edges(ee, 1);
    jj = edges(ee, 2);
    lambda_e = lambdas(ee);
    tij_e = tijs(:, ee);
    T_i = T(:, ii);
    T_j = T(:, jj);
    R_i = R(:, :, ii);
    a = T_i - T_j;
    b = R_i * tij_e;
    cost_ee = trace(a' * a + 2 * lambda_e * (a' * b) + lambda_e^2 * (b' * b)); 

    % if 1-a_log*lambda_e > 0
    %     scale_compensation_ee = -log(1-a_log*lambda_e) - a_log*lambda_e + lambda_e * lambda_e;
    % else
    %     scale_compensation_ee = 0.0;
    % end

    l = lambda_e;
    if l <= 1/a_log
        scale_compensation_ee = 1e+10;
    elseif l<1
        scale_compensation_ee=-1/a_log*log(a_log*l-1)...
            +1/(a_log-1)*(l-1)...
            +b_log/2*(l-1)^2;
    else
        scale_compensation_ee=0;
    end

    cost_out = cost_out + cost_ee + rho * scale_compensation_ee;
end

end %file function