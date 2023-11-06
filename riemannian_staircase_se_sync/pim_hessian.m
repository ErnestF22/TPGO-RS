function v_max = pim_hessian(x, f, thresh)
%PIM_HESSIAN (PIM is an acronym for Power Iteration Method) 
% Iterative method that processes function f at point x in a similar fashion
% to applying the power augmentation method with the linear map f in place
% of the matrix of which the maximum eigenvalue would be computed

if ~exist('thresh','var')
  thresh=1e-5;
end

iterative_change = 1e+6;
d = size(A, 1);
N = size(A, 3);
v = rand(d, 1);



hes = f(x, v);

iteration_num = 0;
while (iterative_change > thresh) && (iteration_num < 1000)

    iteration_num = iteration_num + 1;
    v_prev = v;
    v = v ./ norm(v);
    v = - A * v;
    iterative_change = min(norm(v_prev - v), norm(v_prev + v));
end

v_max = v;

end

