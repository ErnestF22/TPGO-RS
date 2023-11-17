function [lambda_max, v_max] = pim_hessian(x, problem, thresh)
%PIM_HESSIAN (PIM is an acronym for Power Iteration Method) 
% Iterative method that processes function f at point x in a similar fashion
% to applying the power augmentation method with the linear map f in place
% of the matrix of which the maximum eigenvalue would be computed

if ~exist('thresh','var')
  thresh=1e-5;
end

iterative_change = 1e+6;

v_prev = stiefel_randTangentNormVector(x);
v_prev = stiefel_normalize(x, v_prev);
v = v_prev;

iteration_num = 0;
while (iterative_change > thresh) && (iteration_num < 1000)
    iteration_num = iteration_num + 1;
    v_prev = v;
    v = stiefel_normalize(x, v);
    v = - som_rhess_rot_stiefel(x, v_prev, problem);
    iterative_change = min( max(stiefel_metric(x, v_prev, v, 'canonical')), ...
        max(stiefel_metric(x, v, v_prev, 'canonical')) ) ;
end

v_max = stiefel_normalize(x, v);
lambda_max = sum(stiefel_metric([], v_max, som_rhess_rot_stiefel(x, v_max, problem))) / ...
    sum(stiefel_metric([], v_max, v_max));


end

