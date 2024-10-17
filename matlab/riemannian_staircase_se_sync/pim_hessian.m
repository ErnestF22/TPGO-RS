function [lambda_max, v_max] = pim_hessian(x, problem_struct, thresh)
%PIM_HESSIAN (PIM is an acronym for Power Iteration Method) 
% Iterative method that processes function f at point x in a similar fashion
% to applying the power augmentation method with the linear map f in place
% of the matrix of which the maximum eigenvalue would be computed

if ~exist('thresh','var')
  thresh=1e-5;
end

manif = stiefelfactory(problem_struct.sz(1),problem_struct.sz(2),problem_struct.sz(3));
problem_manopt.M = manif;

iterative_change = 1e+6;

% v_start = make_rand_stiefel_3d_array(size(x,1),size(x,2),size(x,3));
% v_start_proj = stiefel_projection_tmp(x, v_start);
% v_start_proj_norm = problem_manopt.M.norm(x, v_start_proj);

v_start_proj_norm = rand_stiefel_vec_tmp(x, size(x,1),size(x,2),size(x,3), "double");
v = v_start_proj_norm;

iteration_num = 0;
while (iterative_change > thresh) && (iteration_num < 1000)
    iteration_num = iteration_num + 1;
    v_prev = v;
    v = v ./ (norm(v(:)));
    v = - som_rhess_rot_stiefel(x, v, problem_struct);
%     iterative_change = min( problem_manopt.M.inner(x, v_prev, v), ...
%         problem_manopt.M.inner(x, v, v_prev) ) ;
end

v_max = v;
lambda_max = problem_manopt.M.inner(x, v_max, som_rhess_rot_stiefel(x, v_max, problem_struct)) / ...
    problem_manopt.M.inner(x, v_max, v_max);


end


function Up = stiefel_projection_tmp(X, U)
        
XtU = multiprod(multitransp(X), U);
symXtU = multisym(XtU);
Up = U - multiprod(X, symXtU);

end

function U = rand_stiefel_vec_tmp(X, n, p, k, array_type)
    U = stiefel_projection_tmp(X, randn(n, p, k, array_type));
    U = stiefel_normalize(X,U);
end

