function [R, R_cost, R_info, R_options] = som_stepone_manopt_stiefel_boumal(T_globalframe, Tijs_vec, edges, fixed_cost_term, num_rows_stiefel, params, R_initguess)
%SOM_STEPONE_MANOPT_STIEFEL Do first step of Manopt pipeline (rotation 
% estimation) on Stiefel manifold

d = params.d;
N = params.N;
initguess_is_available = params.initguess_is_available;


T_globalframe_stiefel = cat_zero_row(T_globalframe, num_rows_stiefel-d);

[L_stiefel, P_stiefel] = make_LT_PT_noloops_stiefel(T_globalframe_stiefel, Tijs_vec, edges, num_rows_stiefel, params);


% R_manopt_stiefel = stiefelstackedfactory(N,d,d);
R_manopt_stiefel = stiefelfactory(num_rows_stiefel,d,N);
fprintf("R_manopt_stiefel size %g\n", R_manopt_stiefel.dim());

step1.M = R_manopt_stiefel; %M = manifold

% Define the step1 cost function and its gradient.
step1.cost = @(x) mycost(x, L_stiefel, P_stiefel, fixed_cost_term);

%step1.egrad = @(x) (L_T+L_T')*x + P_T;
step1.egrad = @(x) myeuclgradient(x, L_stiefel, P_stiefel);
%     step1.grad = @(x) R_manopt_stiefel.egrad2rgrad(matStack(x), step1.egrad);
%     %egrad2rgrad does not seem to work in this case...

% disp(params.riem_grad_mode)
if(strcmp(params.riem_grad_mode, 'manual') == 1)
    step1.grad = @(x) step1.M.proj(x, myeuclgradient(x, L_stiefel, P_stiefel));
    fprintf("Rot grad mode: manual\n");
else
    fprintf("Rot grad mode: auto\n");
end
if(strcmp(params.hessian_mode, 'manual') == 1)
    step1.hess = @(x, u) step1.M.proj(x, myeuclhess(x, u, L_stiefel, P_stiefel));
    fprintf("Rot Hess mode: manual\n");
else
    fprintf("Rot Hess mode: auto\n");
end
% Numerically check gradient consistency (optional).
%checkgradient(step1);

% Solve
% options.maxiter = 5;
options.verbosity = 0;
if initguess_is_available
    disp("Manopt_stiefel eval rot with initguess");
    [R, R_cost, R_info, R_options] = trustregions(step1, R_initguess, options);
else
    disp("Manopt_stiefel eval rot WITHOUT initguess");
    [R, R_cost, R_info, R_options] = trustregions(step1, [], options);
end

end %file function

%%
function f = mycost(x, L_stiefel, P_stiefel, fixed_cost_term)
f = trace(matStack(x)' * L_stiefel * matStack(x) + matStack(x)' * P_stiefel) + fixed_cost_term;
end
%%
function g = myeuclgradient(x, L_stiefel, P_stiefel)
g = matUnstack(L_stiefel*matStack(x) + (L_stiefel')*matStack(x) + P_stiefel, size(x, 1));
end
%%
% function h = myriemgradient(x, L_stiefel, P_stiefel)
% g = myeuclgradient(x,L_stiefel,P_stiefel);    
% h = multiprod(multitransp(x), g) - multiprod(multitransp(g), x);
% h = 0.5 .* h;
% end
%%
function eucl_hess = myeuclhess(x, u, L, P)
P_3d = matUnstack(P, size(x, 1));
eucl_hess = multiprod3(u, multitransp(x), P_3d) ... 
        + multiprod3(x, multitransp(u), P_3d) ...
        - multiprod3(u, multitransp(P_3d), x) ...
        - multiprod3(x, multitransp(P_3d), u);
eucl_hess = 0.5 .* eucl_hess;
end
%%
% function hess = myhess(x,u,L,P)
% g = myeuclhess(x,u,L,P);
% hess = multiprod(multitransp(x), g) - multiprod(multitransp(g), x);
% hess = 0.5.*hess;
% end
%%
