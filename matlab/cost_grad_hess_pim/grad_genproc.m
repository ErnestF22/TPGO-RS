% Riemannian gradient of the cost function.
function [g] = grad_genproc(X, problem_data)
    R = X.R;
    T = X.T;
    Tijs = problem_data.Tijs;
    edges = problem_data.edges;
    [problem_structs.P, problem_structs.frct] = ...
        make_step1_p_fct(T, Tijs, edges);
    [problem_structs.LR, problem_structs.PR, problem_structs.BR] = ...
        make_LR_PR_BR_noloops(R, Tijs, edges);
    g.R = rsom_rgrad_rot_stiefel(R, problem_structs);
    g.T = rsom_rgrad_transl_stiefel(T, problem_structs);
end %grad