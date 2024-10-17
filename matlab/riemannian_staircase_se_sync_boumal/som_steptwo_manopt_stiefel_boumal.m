function [T, T_cost, T_info, T_options] = som_steptwo_manopt_stiefel_boumal(R_stiefel, T_globalframe, Tijs_vec, edges, num_rows_stiefel, params, T_initguess)
%SOM_STEPTWO_MANOPT Do second step of Manopt pipeline (translation estimation)
    
    N = params.N;
    d = params.d;
    initguess_is_available = params.initguess_is_available;

    T_globalframe_stiefel = cat_zero_row(T_globalframe, num_rows_stiefel-d);

    %     R_stiefel = zeros(num_rows_stiefel, d, N);
    %     for ii = 1:size(R, 3)
    %         R_stiefel(:,:,ii) = cat_zero_row(R(:,:,ii));
    %     end
    %     R_2d_stiefel = matStack(R_truth_stiefel);

    [A,b] = make_A_b_noloops_stiefel(R_stiefel, T_globalframe_stiefel, Tijs_vec, edges, num_rows_stiefel, params);

    % R_manopt = stiefelstackedfactory(N,d,d);
    transl_manopt = euclideanfactory(num_rows_stiefel*N,1);
    fprintf("R_manopt size %g\n", transl_manopt.dim());
    
    step1.M = transl_manopt; %M = manifold
    
    % Define the step1 cost function and its Euclidean gradient.
    step1.cost = @(x) mycost(x, A, b);

    % !! since the egrad2rgrad comparison is not necessary for the
    % euclidean factory, the if('manual') is better to be avoided
%     if(strcmp(params.hessian_mode, 'manual') == 1)
    step1.grad = @(x) myeuclgradient(x, A, b);
%         fprintf("Transl grad mode: manual\n");
%     else
%         fprintf("Transl grad mode: auto\n");
%     end
    
    if(strcmp(params.hessian_mode, 'manual') == 1)
        step1.hess = @(x, u) myhess(x, u, A, b);
        fprintf("Transl Hess mode: manual\n");
    else
        fprintf("Transl Hess mode: auto\n");
    end

    % Numerically check gradient consistency (optional).
    % checkgradient(step1);
    
    %     options.maxiter = 5;
    options.verbosity = 0;

    % Solve.
    if initguess_is_available
        disp("Manopt eval transl with initguess");
        [T, T_cost, T_info, T_options] = trustregions(step1, T_initguess(:), options);
    else
        disp("Manopt eval transl WITHOUT initguess");
        [T, T_cost, T_info, T_options] = trustregions(step1, [], options);
    end
end


%%
function f = mycost(x, A_stiefel, b)
    e = A_stiefel*x(:) + b;
    f = e'*e;
end
%%
function g = myeuclgradient(x, A_stiefel, b) 
    g = 2 * (A_stiefel' * A_stiefel) * x(:) + 2*A_stiefel'*b;
end
%%
% function h = myriemgradient(x, L, P) 
%     %g = myeuclgradient(x,L,P) BUT still stacked
%     g = L*matStack(x) + (L')*matStack(x) + P;
%     h = zeros(size(g));
%     d = size(g,2);
%     for ii = 1:d:size(g,1)
%         gpart = g(ii:ii+d-1, :);
%         g_asym = 0.5 .* (gpart - gpart');
%         h(ii:ii+d-1, :) = g_asym' * gpart';
%     end
%     h = matUnstack(h);
% end
%%
function eucl_hess = myhess(x,u,A_stiefel,b)
eucl_hess = 2.*(A_stiefel'*A_stiefel)*u;
end

