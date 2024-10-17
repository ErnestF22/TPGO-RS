function [lambda_max, x_max] = pim_function_genproc(f, x_start, normalization_fun, thresh)
%PIM_FUNCTION (where PIM is an acronym for Power Iteration Method) 
% Iterative method that returns an eigenvector associated to the maximum 
% eigenvalue of matrix A
% x_start is the random vector from which PIM starts

if ~exist('thresh','var')
  thresh=1e-10;
end

%% R
iterative_change = 1e+6;
xR = x_start.R;
xT = x_start.T;
xfull = [matStackH(x_start.R), x_start.T];
iteration_num = 0;
while (iteration_num < 2500) % && (abs(iterative_change) > thresh)
    iteration_num = iteration_num + 1;
    x_prev_R = xR;
    x_prev_T = xT;
    xfull_prev = [matStackH(x_prev_R), x_prev_T];
    norm_RT = norm(xfull);
    x.R = xR/norm_RT;
    x.T = xT/norm_RT;
    fx = f(x);
    xR = -fx.R;
    xT = -fx.T;
    xfull = [matStackH(xR), xT];
    iterative_change = max(normalization_fun(xfull_prev - xfull), [],"all");
end

norm_RT_max = norm([matStackH(x.R), x.T]);
x_max.R = x.R/norm_RT_max;
x_max.T = x.T/norm_RT_max;
f_x_max = f(x_max);
% lambda_max_R = sum(stiefel_metric([], (x_max.R), f_x_max.R)) / ...
%     sum(stiefel_metric([], x_max.R, x_max.R));
% lambda_max_T = sum(stiefel_metric([], (x_max.T), f_x_max.T)) / ...
%     sum(stiefel_metric([], x_max.T, x_max.T));
lambda_max = x_max.R(:)' * f_x_max.R(:) + x_max.T(:)' * f_x_max.T(:);
full_xmax = [matStackH(x_max.R), x_max.T];
lambda_max = lambda_max / ...
    sum(stiefel_metric([], full_xmax(:), full_xmax(:))); %note: simple eucl. metric

end %file function

