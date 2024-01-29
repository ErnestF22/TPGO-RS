function [R_out, R_cost, R_info, R_options] = rsom_step1( ...
    T_gf, Tijs, edges, params, R_initguess)
%DO_RSOM_MANOPT_STEP1 Run step 1 of the Manopt pipeline and return Manopt
% optimization outputs.



%problem size params
nrs = size(T_gf, 1);
d = size(Tijs, 1);
N = size(T_gf, 2);
problem_step1.sz = [nrs, d, N];
% num_edges = size(edges, 1);

if nargin < 5
    R_initguess = make_rand_stiefel_3d_array(nrs, d, N);
end

% Create the problem structure.
manifold = stiefelfactory(nrs,d,N);
problem_step1.M = manifold;

[P, frct] = make_step1_p_fct(T_gf, Tijs, edges);
problem_step1.P = P;
problem_step1.frct = frct;
 
% Define the problem cost function and its Euclidean gradient.
problem_step1.cost  = @(x) rsom_cost_rot_stiefel(x,problem_step1);
problem_step1.egrad = @(x) rsom_egrad_rot_stiefel(x,problem_step1);
problem_step1.grad = @(x) rsom_rgrad_rot_stiefel(x,problem_step1);
% problem.ehess = @(x,u) rsom_ehess_rot_stiefel(x,u,problem);
% problem.hess = @(x,u) rsom_rhess_rot_stiefel(x,problem);
 
% Numerically check gradient and hessian consistency (optional).
% checkgradient(problem);
% checkhessian(problem);
 
% Solve.
if params.initguess_is_available
    [R_out, R_cost, R_info, R_options] = trustregions(problem_step1, R_initguess);
else
    [R_out, R_cost, R_info, R_options] = trustregions(problem_step1);
end


end

