%%
% It is not necessary to define the Hessian of the cost. We do it
% mostly to illustrate how to do it and to study the spectrum of the
% Hessian at the solution (see further down).
function h_sh = ssom_rhess_genproc_shifted(X, Xdot, mu_shift, problem_data)
    % R = X.R;
    % T = X.T;
    % lambdas = X.lambda;
    Rdot = Xdot.R;
    Tdot = Xdot.T;
    lambdasdot = Xdot.lambda;

    h = ssom_rhess_genproc(X, Xdot, problem_data);

    h_sh.R = h.R - mu_shift*Rdot;
    h_sh.T = h.T - mu_shift*Tdot;
    h_sh.lambda = h.lambda - mu_shift*lambdasdot;
end %hess
