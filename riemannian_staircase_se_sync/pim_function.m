function [lambda_max, x_max] = pim_function(f, x_start, normalization_fun, thresh)
%PIM_FUNCTION (where PIM is an acronym for Power Iteration Method) 
% Iterative method that returns an eigenvector associated to the maximum 
% eigenvalue of matrix A
% x_start is the random vector from which PIM starts

if ~exist('thresh','var')
  thresh=1e-10;
end

iterative_change = 1e+6;
dim = size(x_start);
x = x_start;
iteration_num = 0;
iterative_change = 1e+6;
while (iteration_num < 3000) % && (abs(iterative_change) > thresh)
    iteration_num = iteration_num + 1;
    x_prev = x;
    x = normalization_fun(x);
    x = - f(x);
    iterative_change = max(min(normalization_fun(x_prev - x), ...
        normalization_fun(x_prev + x)),[],"all");
end

x_max = normalization_fun(x);
lambda_max = sum(stiefel_metric([], (x_max), f(x_max))) / ...
    sum(stiefel_metric([], x_max, x_max));

end

