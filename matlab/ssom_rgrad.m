function g = ssom_rgrad(X, problem_data)
    %g.R
    g.R = rgrad_R(X, problem_data);
    %g.T
    g.T = egrad_T(X, problem_data);
    %g.lambda
    g.lambda = ssom_grad_lambda(X, problem_data);

end

function g=rgrad_R(X, problem_data)

tijs_scaled = make_tijs_scaled(X.lambda, problem_data.tijs);
[P, ~] = make_step1_p_fct(X.T, tijs_scaled, problem_data.edges);

R = X.R;
d = size(R, 2);
eg=matUnstackH(P,d);

g = stiefel_tangentProj(R, eg);
end

function g=egrad_T(X,problem_data)
T = X.T;
tijs_scaled = make_tijs_scaled(X.lambda, problem_data.tijs);
[LR, PR, ~] = make_LR_PR_BR_noloops(X.R, tijs_scaled, problem_data.edges);
g=T*(LR+LR')+(PR)';
end

