function Y_out = linesearch_decrease_hessian(problem, x, v, initial_cost, alphas)
% LINESEARCH_DECREASE_HESSIAN
% Backtracking line-search aiming merely for a decrease in cost value.
%
% The Line-search algorithm based on a simple backtracking method. The
% search direction v provided has to be a descent direction, but needs not to
% be a first-order descent,
% i.e.: this line-search can be used even if x is a
% critical point, as long as the cost function is strictly decreasing
% along the direction d.
%
% The line-search merely guarantees a decrease in the cost (unless a
% stopping criterion triggers first, such as exceeding a maximal number of
% iterations). This is typically useful to escape saddle points (critical
% points admitting descent directions at the second order). Escape
% directions can be computed using the Hessian eigenvectors corresponding
% to negative eigenvalues, for example.
%


if ~exist('alphas','var')
    alphas=linspace(-1e-3, 1e-3, 21); %t
end



retractions_vec = zeros(size(alphas));
plot_vals_taylor = zeros(size(alphas));
found_lower = boolean(0);

pos_start = ceil(length(alphas)/2) + 1; %TODO: check if this works with even alphas size
neg_end = floor(length(alphas)/2) - 1; %TODO: check if this works with even alphas size

% First, go towards 0+ (positive direction)
for ii = pos_start:length(alphas)
    candidate_val = x + alphas(ii) * v;
    %check
    disp("Is candidate_val on Stiefel?")
    if check_is_on_stiefel(candidate_val)
        disp("YES")
    else
        disp("NO");
    end
    cand_val_cost = problem.cost(candidate_val);
    disp("[initial_cost, cand_val_cost]")
    disp([initial_cost, cand_val_cost])


    plot_vals_taylor(ii) = initial_cost+...
        alphas(ii)^2/2*sum(stiefel_metric(x,v,problem.hess(x,v),'euclidean'));
    if plot_vals_taylor(ii) < initial_cost
        found_lower = boolean(1);
        Y_out = candidate_val;
        break;
    end
end

if ~found_lower
    % Else, go towards 0+ (positive direction)
    for ii = 1:neg_end
        candidate_val = x + alphas(ii) * v;
        %check
        disp("Is candidate_val on Stiefel?")
        if check_is_on_stiefel(candidate_val)
            disp("YES")
        else
            disp("NO");
        end
        plot_vals_taylor(ii) = initial_cost+...
            alphas(ii)^2/2*sum(stiefel_metric(x,v,problem.hess(x,v),'euclidean'));
        if plot_vals_taylor(ii) < initial_cost
            found_lower = boolean(1);
            Y_out = candidate_val;
            break;
        end
    end
end

%linesearch decrease has failed
Y_out = x;



end %file function