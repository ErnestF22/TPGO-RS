% Define the cost function here. Points on the manifold M are
% structures with fields X.A and X.R, containing matrices of sizes
% respectively nxm and nxnxN. 
% The store structure (the caching system) is used to keep the residue 
% matrix E in memory, as it is also used in the computation of the 
% gradient and of the Hessian. This way, we prevent redundant computations.
function [f] = cost_genproc(X, problem_data)
%     f = rsom_cost_rot_stiefel(R, problem);
    f = rsom_cost_base(X, problem_data);
end %cost   