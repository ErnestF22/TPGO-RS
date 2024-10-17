function [R, R_cost, R_info, R_options] = som_stepone_manopt(T_globalframe, Tijs_vec, edges, fixed_cost_term, params, R_initguess)
%SOM_STEPONE_MANOPT Do first step of Manopt pipeline (rotation estimation)

d = params.d;
N = params.N;
initguess_is_available = params.initguess_is_available;

[L, P] = make_LT_PT_noloops(T_globalframe, Tijs_vec, edges, params);

% R_manopt = stiefelstackedfactory(N,d,d);
R_manopt = rotationsfactory(d,N);
fprintf("R_manopt size %g\n", R_manopt.dim());

step1.M = R_manopt; %M = manifold

% Define the step1 cost function and its Euclidean gradient.
step1.cost = @(x) mycost(x, L, P, fixed_cost_term);

step1.egrad = @(x) myeuclgradient(x, L, P);

% disp(params.riem_grad_mode)
if(strcmp(params.riem_grad_mode, 'manual') == 1)
    step1.grad = @(x) myriemgradient(x, L, P);
    fprintf("Rot grad mode: manual\n");
else
    fprintf("Rot grad mode: auto\n");
end
if(strcmp(params.hessian_mode, 'manual') == 1)
    step1.hess = @(x, u) myhess(x, u, L, P);
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
    disp("Manopt eval rot with initguess");
    [R, R_cost, R_info, R_options] = trustregions(step1, R_initguess, options);
else
    disp("Manopt eval rot WITHOUT initguess");
    [R, R_cost, R_info, R_options] = trustregions(step1, [], options);
end

end %file function

%%
function f = mycost(x, L, P, fixed_cost_term)
f = trace((matStack(x))'*L*(matStack(x)) + (matStack(x))'*P) + fixed_cost_term;
end
%%
function g = myeuclgradient(x, L, P)
g = matUnstack(L*matStack(x) + (L')*matStack(x) + P);
end
%%
function h = myriemgradient(x, L, P)
g = myeuclgradient(x,L,P);    
h = multiprod(multitransp(x), g) - multiprod(multitransp(g), x);
h = 0.5 .* h;
end
%%
function eucl_hess = myeuclhess(x, u, L, P)
P_3d = matUnstack(P);
eucl_hess = multiprod3(u, multitransp(x), P_3d) ... 
        + multiprod3(x, multitransp(u), P_3d) ...
        - multiprod3(u, multitransp(P_3d), x) ...
        - multiprod3(x, multitransp(P_3d), u);
eucl_hess = 0.5 .* eucl_hess;
end
%%
function hess = myhess(x,u,L,P)
g = myeuclhess(x,u,L,P);
hess = multiprod(multitransp(x), g) - multiprod(multitransp(g), x);
hess = 0.5.*hess;
end
%%