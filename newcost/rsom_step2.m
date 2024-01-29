function [T_out, T_cost, T_info, T_options] = rsom_step2( ...
    R_gf, Tijs, edges, params)
%DO_RSOM_MANOPT_STEP1 Run step 1 of the Manopt pipeline and return Manopt
% optimization outputs.

%problem size params
nrs = size(R_gf, 1);
d = size(R_gf, 2);
N = size(R_gf, 3);
problem_step2.sz = [nrs, d, N];
% num_edges = size(edges, 1);

% Create the problem structure.
manifold = euclideanfactory(nrs,N);
problem_step2.M = manifold;

[LR,PR,BR] = make_LR_PR_BR_noloops( ...
    R_gf, Tijs, edges);
problem_step2.LR = LR;
problem_step2.PR = PR;
problem_step2.BR = BR;
 
% Define the problem cost function and its Euclidean gradient.
problem_step2.cost  = @(x) rsom_cost_transl_stiefel(x,problem_step2);
problem_step2.egrad = @(x) rsom_egrad_transl_stiefel(x,problem_step2);
% problem_step2.grad = @(x) rsom_rgrad_transl_stiefel(x,problem_step2);
problem_step2.ehess = @(x,u) rsom_ehess_transl_stiefel(x,u,problem_step2);
% problem.hess = @(x,u) rsom_rhess_transl_stiefel(x,problem_step2);
 
% Numerically check gradient and hessian consistency (optional).
% checkgradient(problem);
% checkhessian(problem);
 
% Solve.
[T_out, T_cost, T_info, T_options] = trustregions(problem_step2);


end

