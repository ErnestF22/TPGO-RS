function stillzero = check_yplus_grad(Yopt, problem_data, measurements, num_rows_stiefel, som_params, eps)
%FUNCTION CHECK_YPLUS_GRAD() Returns boolean stillzero that states whether
%max(Yplus, [], 'all') is close enough to zero (i.e., < eps)

if (nargin < 6)
    eps = 1e-5;
end

R_manopt_stiefel = stiefelfactory(num_rows_stiefel,problem_data.d,problem_data.n);
fprintf("R_manopt_stiefel size %g\n", R_manopt_stiefel.dim());
manopt_data.M = R_manopt_stiefel; %M = manifold
% myeuclgradient(Yopt, problem_data, measurements, som_params)
% manopt_data.grad = manopt_data.M.proj(Yopt, myeuclgradient(Yopt, problem_data, measurements, som_params))
Yopt_grad = manopt_data.M.proj(Yopt, myeuclgradient(Yopt, problem_data, measurements, som_params));

Yplus = cat_zero_rows_3d_array(Yopt);
R_manopt_stiefel = stiefelfactory(num_rows_stiefel+1,problem_data.d,problem_data.n);
fprintf("R_manopt_stiefel size %g\n", R_manopt_stiefel.dim());
manopt_data.M = R_manopt_stiefel; %M = manifold
% myeuclgradient(Yplus, problem_data, measurements, som_params)
% manopt_data.grad = manopt_data.M.proj(Yplus, myeuclgradient(Yplus, problem_data, measurements, som_params))
Yplus_grad = manopt_data.M.proj(Yplus, myeuclgradient(Yplus, problem_data, measurements, som_params));

stillzero = max(abs(Yplus_grad), [], 'all') < eps;

if (stillzero == false)
   disp("Yopt_grad");
   disp(Yopt_grad);
   disp("Yplus_grad");
   disp(Yplus_grad);
end

end %file function


%%
function f = mycost(x, problem_data, measurements, som_params)
    measurements.T_globalframe_stiefel = cat_zero_row(measurements.T_globalframe_stiefel, measurements.d_stiefel-problem_data.d);
    [L_stiefel, P_stiefel] = make_LT_PT_noloops_stiefel( ...
        measurements.T_globalframe_stiefel, measurements.Tijs_vec, measurements.edges, measurements.d_stiefel, som_params);
    cost_const_term_tij = compute_fixed_cost_term(measurements.Tijs_vec, problem_data.d);
    f = trace(matStack(x)' * L_stiefel * matStack(x) + matStack(x)' * P_stiefel) + cost_const_term_tij;
end
%%
function g = myeuclgradient(x, problem_data, measurements, som_params)
    measurements.T_globalframe_stiefel = cat_zero_row(measurements.T_globalframe_stiefel, size(x,1)-problem_data.d);
    [L_stiefel, P_stiefel] = make_LT_PT_noloops_stiefel( ...
        measurements.T_globalframe_stiefel, measurements.Tijs_vec, measurements.edges, size(x,1), som_params);
    g = matUnstack(L_stiefel*matStack(x) + (L_stiefel')*matStack(x) + P_stiefel, size(x, 1));
end
%%
% function h = myriemgradient(x, L_stiefel, P_stiefel)
% g = myeuclgradient(x,L_stiefel,P_stiefel);
% h = multiprod(multitransp(x), g) - multiprod(multitransp(g), x);
% h = 0.5 .* h;
% end
%%
function eucl_hess = myeuclhess(x, u, problem_data, measurements, som_params)
measurements.T_globalframe_stiefel = cat_zero_row(measurements.T_globalframe_stiefel, size(x,1)-problem_data.d);
[L_stiefel, P_stiefel] = make_LT_PT_noloops_stiefel( ...
        measurements.T_globalframe_stiefel, measurements.Tijs_vec, measurements.edges, size(x,1), som_params);
P_3d = matUnstack(P_stiefel, size(x, 1));
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
