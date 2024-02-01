function [Y0, lambda_pim, v_pim] = rsom_pim_hessian(R, problem_struct_next, thresh)
%RSOM_PIM_HESSIAN Summary of this function goes here
%   Detailed explanation goes here

%!!! Only with NON-shifting PIM

rhess_fun_han = @(u) rsom_rhess_rot_stiefel(R,u,problem_struct_next);

stiefel_normalize_han = @(x) x./ (norm(x(:)));

u_start = stiefel_randTangentNormVector(R);
[lambda_pim, v_pim] = pim_function(rhess_fun_han, u_start, stiefel_normalize_han, thresh);
disp('Difference between lambda*v_max and H(v_max) should be in the order of the tolerance:')
eigencheck_hessian(lambda_pim, v_pim, rhess_fun_han);

[~, Y0] = linesearch_decrease(problem_struct_next, ...
    R, -v_pim, rsom_cost_rot_stiefel(R,problem_struct_next));


end %file function

