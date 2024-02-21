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
x = x_start;
iteration_num = 0;
while (iteration_num < 1200) % && (abs(iterative_change) > thresh)
%     x_R = x.R(1:5);
    iteration_num = iteration_num + 1;
    x_prev = x;
    x.R = normalization_fun(x.R);
    x.T = normalization_fun(x.T);
    fx = f(x);
    x.R = -fx.R;
    x.T = -fx.T;
    iterative_change = max(min(normalization_fun(x_prev.R - x.R), ...
        normalization_fun(x_prev.R + x.R)),[],"all");
end

x_max.R = normalization_fun(x.R);
x_max.T = normalization_fun(x.T); %!! used in f(x_max)
f_x_max = f(x_max);
lambda_max.R = sum(stiefel_metric([], (x_max.R), f_x_max.R)) / ...
    sum(stiefel_metric([], x_max.R, x_max.R));
lambda_max.T = sum(stiefel_metric([], (x_max.T), f_x_max.T)) / ...
    sum(stiefel_metric([], x_max.T, x_max.T));

% x_max_cand_R = x_max.R;
% lambda_max_cand_R = lambda_max.R;
% 
% %% T
% iterative_change = 1e+6;
% x = x_start;
% iteration_num = 0;
% while (iteration_num < 200) % && (abs(iterative_change) > thresh)
%     iteration_num = iteration_num + 1;
%     x_prev = x;
%     x.T = normalization_fun(x.T);
%     fx = f(x);
%     x.R = -fx.R;
%     x.T = -fx.T;
%     iterative_change = max(min(normalization_fun(x_prev.T - x.T), ...
%         normalization_fun(x_prev.T + x.T)),[],"all");
% end
% 
% x_max.R = normalization_fun(x.R); %!! used in f(x_max)
% x_max.T = normalization_fun(x.T);
% f_x_max = f(x_max);
% lambda_max.T = sum(stiefel_metric([], (x_max.T), f_x_max.T)) / ...
%     sum(stiefel_metric([], x_max.T, x_max.T));
% 
% %%
% x_max.R = x_max_cand_R;
% lambda_max.R = lambda_max_cand_R;
% % x_max.T = normalization_fun(x);
% % lambda_max.T = sum(stiefel_metric([], (x_max.T), f(x_max.T))) / ...
% %     sum(stiefel_metric([], x_max, x_max));

end %file function

