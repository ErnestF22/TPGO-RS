%%
% It is not necessary to define the Hessian of the cost. We do it
% mostly to illustrate how to do it and to study the spectrum of the
% Hessian at the solution (see further down).
function hLambda = ssom_hess_genproc_lambda(X, Xdot, problem_data)
    h = ssom_rhess_genproc(X, Xdot, problem_data);
    hLambda = h.lambda;
end %hess
