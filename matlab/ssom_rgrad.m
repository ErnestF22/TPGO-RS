function g = ssom_rgrad(x, problem_data)
    lambdas = x.lambda;
    T = x.T;
    R = x.R;

    edges = problem_data.edges;

    tijs_scaled = make_tijs_scaled(lambdas, problem_data.Tijs);
    
    %g.R
    problem_data_R = problem_data;
    problem_data_R.Tijs = tijs_scaled;
    [problem_data_R.P, problem_data_R.frct] = ...
        make_step1_p_fct(T, tijs_scaled, edges);
    g.R = rsom_rgrad_rot_stiefel(R, problem_data_R);
    %g.T
    problem_data_T = problem_data;
    problem_data_T.Tijs = tijs_scaled;
    [problem_data_T.LR, problem_data_T.PR, problem_data_T.BR] = ...
        make_LR_PR_BR_noloops(R, tijs_scaled, edges);
    g.T = rsom_rgrad_transl_stiefel(T, problem_data_T);
    %g.lambda
    problem_data_lambdas = problem_data;
    problem_data_lambdas.T = T;
    problem_data_lambdas.R = R;
    g.lambda = ssom_grad_lambda(lambdas, problem_data_lambdas);

end