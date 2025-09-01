function g = ssom_egrad(X, problem_data)

    R = X.R;
    T = X.T;
    lambdas = X.lambda;
    %g.R
    g.R = ssom_egrad_R(R,T,lambdas, problem_data);
    %g.T
    g.T = ssom_egrad_T(R,T,lambdas, problem_data);
    %g.lambda
    g.lambda = ssom_egrad_lambda(R,T,lambdas, problem_data);

end


