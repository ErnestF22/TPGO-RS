function g = lsom_rgrad_relu(X, problem_data)
    R = X.R;
    T = X.T;
    lambdas = X.lambda;
    % mu = problem_data.mu;
    % y = problem_data.y;
    % z = problem_data.z;
    %g.R
    g.R = ssom_rgrad_R(R, T, lambdas, problem_data); %Note: using this notation as it facilitates geodesic tests
    %g.T
    g.T = ssom_egrad_T(R, T, lambdas, problem_data);
    %g.lambda
    g.lambda = ssom_rgrad_lambda_relu(R, T, lambdas, problem_data);

end