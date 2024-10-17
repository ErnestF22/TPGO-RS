%%
% It is not necessary to define the Hessian of the cost. We do it
% mostly to illustrate how to do it and to study the spectrum of the
% Hessian at the solution (see further down).
function [h] = hess_genproc(X, Xdot, problem_data)
    R = X.R;
    T = X.T;
    Rdot = Xdot.R;
    Tdot = Xdot.T;
    % Careful: tangent vectors on the rotation group are represented as
    % skew symmetric matrices. To obtain the corresponding vectors in
    % the ambient space, we need a little transformation. This
    % transformation is typically not needed when we compute the
    % formulas for the gradient and the Hessian directly in Riemannian
    % form instead of resorting the egrad2rgrad and ehess2rhess. These
    % latter tools are convenient for prototyping but are not always
    % the most efficient form to execute the computations.
    Tijs = problem_data.Tijs;
    edges = problem_data.edges;
    [problem_structs.P, problem_structs.frct] = ...
        make_step1_p_fct(T, Tijs, edges);
    [problem_structs.LR, problem_structs.PR, problem_structs.BR] = ...
        make_LR_PR_BR_noloops(R, Tijs, edges);
 
    hrt = compute_hrt(X,Xdot,Tijs,edges);
    htr = compute_htr(X,Xdot,Tijs,edges);

    h.R = rsom_rhess_rot_stiefel(R, Rdot, problem_structs) + hrt;
    h.T = rsom_rhess_transl_stiefel(T, Tdot, problem_structs) + htr;
end %hess
