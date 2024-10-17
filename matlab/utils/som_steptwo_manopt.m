function [T, T_cost, T_info, T_options] = som_steptwo_manopt(R, T_globalframe_nois, Tijs_vec_nois, edges, params, T_initguess)
%SOM_STEPTWO_MANOPT Do second step of Manopt pipeline (translation estimation)
    
    N = params.N;
    d = params.d;
    initguess_is_available = params.initguess_is_available;
  
    [A,b] = make_A_b_noloops(R, T_globalframe_nois, Tijs_vec_nois, edges, params);

    % R_manopt = stiefelstackedfactory(N,d,d);
    transl_manopt = euclideanfactory(d*N, 1);
    fprintf("R_manopt size %g\n", transl_manopt.dim());
    
    step1.M = transl_manopt; %M = manifold
    
    % Define the step1 cost function and its Euclidean gradient.
    step1.cost = @(x) mycost(x, A, b);
    
    step1.egrad = @(x) myeuclgradient(x, A, b);
%     step1.grad = @(x) R_manopt.egrad2rgrad(matStack(x), step1.egrad);

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
        [T, T_cost, T_info, T_options] = trustregions(step1, T_initguess, options);
    else
        disp("Manopt eval transl WITHOUT initguess");
        [T, T_cost, T_info, T_options] = trustregions(step1, [], options);
    end
end


%%
function f = mycost(x, A, b)
    e=A*x + b;
    f = e'*e;
end
%%
function g = myeuclgradient(x, A, b) 
    g = 2 * (A' * A) * x + 2*A'*b;
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
function eucl_hess = myhess(x,u,A,b)
eucl_hess = 2.*(A'*A)*u;
end
%%
