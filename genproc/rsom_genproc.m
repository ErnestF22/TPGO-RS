function transf_out = rsom_genproc(T_gf, Tijs, edges, params, transf_initguess)
%RSOM_RS Rsom Manopt pipeline, with the addition of the Riemannian
%Staircase ("RS")

if ~exist('thresh','var')
  thr=1e-5;
end

nrs = size(T_gf, 1);
d = size(Tijs, 1);
N = size(T_gf, 2);
problem.sz = [nrs, d, N];

r0 = d+1;
problem_data.Tijs = Tijs;
problem_data.edges = edges;

tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
M = productmanifold(tuple);

% Setup the problem structure with manifold M and cost+grad functions.
problem.M = M;
problem_data.sz = [nrs, d, N];
problem.cost = @(x) cost_genproc(x, problem_data);
problem.grad = @(x) grad_genproc(x, problem_data);
problem.hess = @(x, u) hess_genproc(x, u, problem_data);

% checkgradient(problem);
% checkhessian(problem);

X = trustregions(problem);
T = X.T;
R = X.R;



for staircase_step_idx = r0:d*N+1
    problem_struct_next.sz = [staircase_step_idx, d, N];
    problem_struct_next.Tijs = Tijs;
    problem_struct_next.edges = edges;

    [Y_star, lambda, v] = rsom_pim_hessian_genproc( ...
        X, problem_struct_next, thr);
    disp("v") %just to remove unused variable warning
    disp(v)

    if lambda > 0
        disp("R, T eigenvals > 0: exiting staircase")
        break;
    end
    
    % next optimization iteration
    tuple_next.R = stiefelfactory(staircase_step_idx, d, N);
    tuple_next.T = euclideanfactory(staircase_step_idx, N);
    M_next = productmanifold(tuple_next);
    problem_next.M = M_next;    
    problem_next.cost = @(x) cost_genproc(x, problem_data); %!! problem_data is the same
    problem_next.grad = @(x) grad_genproc(x, problem_data);
    problem_next.hess = @(x, u) hess_genproc(x, u, problem_data);

    
    X = trustregions(problem_next, Y_star);
    T = X.T;
    R = X.R;

    if rank(matStackH(Y_star.R))<staircase_step_idx
        break;
    end

end

%TODO: Improve reprojection on SE(d)^N

% METHOD 1): Simply extract upper left submats
% transf_out = RT2G(R(1:d, 1:d, 1:N), T(1:d, 1:N));

% METHOD 2): SE-SYNC round solution (only for rotations!)
% R_out = matUnstackH( ...
%     round_solution_se_sync(matStackH(R), problem_struct_next));
% T_out = T(1:d, 1:N); %what to do??
% transf_out = RT2G(R_out, T_out);

% METHOD 3): CODEMETA lowRankLocalization_solution_extractProjection
[R_out, T_out] = lowRankLocalization_solution_extractProjection( ...
    matStack(multitransp(R)) * T);
transf_out = RT2G(R_out, T_out);

end %file function

