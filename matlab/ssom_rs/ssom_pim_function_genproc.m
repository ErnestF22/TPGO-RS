function [lambda_max, x_max] = ssom_pim_function_genproc(f, x_start, normalization_fun, thresh)
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
xLambda = x_start.lambda;
xfull = [x_start.R(:); x_start.T(:); x_start.lambda(:)];
iteration_num = 0;
while (iteration_num < 2000) % && (abs(iterative_change) > thresh)
    iteration_num = iteration_num + 1;
    x_prev_R = xR;
    x_prev_T = xT;
    x_prev_lambda = xLambda;
    xfull_prev = [x_prev_R(:); x_prev_T(:); x_prev_lambda(:)];
    norm_RTlambda = norm(xfull);
    x.R = xR/norm_RTlambda;
    x.T = xT/norm_RTlambda;
    x.lambda = xLambda / norm_RTlambda;
    fx = f(x);
    xR = -fx.R;
    xT = -fx.T;
    xLambda = -fx.lambda;
    xfull = [xR(:); xT(:); xLambda(:)];
    % iterative_change = max(normalization_fun(xfull_prev - xfull), [], "all");
end

norm_RTlambda_max = norm([x.R(:); x.T(:); x.lambda(:)]);
x_max.R = x.R/norm_RTlambda_max;
x_max.T = x.T/norm_RTlambda_max;
x_max.lambda = x.lambda/norm_RTlambda_max;
f_x_max = f(x_max);
% lambda_max_R = sum(stiefel_metric([], (x_max.R), f_x_max.R)) / ...
%     sum(stiefel_metric([], x_max.R, x_max.R));
% lambda_max_T = sum(stiefel_metric([], (x_max.T), f_x_max.T)) / ...
%     sum(stiefel_metric([], x_max.T, x_max.T));
lambda_max = x_max.R(:)' * f_x_max.R(:) + ...
    x_max.T(:)' * f_x_max.T(:) + ...
    x_max.lambda(:)' * f_x_max.lambda(:);
full_xmax = [x_max.R(:); x_max.T(:); x_max.lambda(:)];
lambda_max = lambda_max / ...
    sum(stiefel_metric([], full_xmax(:), full_xmax(:))); %note: simple eucl. metric

end %file function

