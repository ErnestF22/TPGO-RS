close all
clear
clc

resetRands(2);

d = 3;
nrs = 4;
N = 5;

% %graph random init
% e = 29;
% G = graph(true(N), 'omitselfloops'); % Alternative without self-loops
% p = randperm(numedges(G), e);
% G = graph(G.Edges(p, :));
% edges = table2array(G.Edges);
% num_edges = e;

%testdata random init
testdata = testNetwork_som(3); %4 is the default option
edges = (testdata.E);
num_edges = testdata.NEdges;
problem_data.edges = edges;
% problem_data.num_edges = num_edges;

R_gf = make_rand_stiefel_3d_array(nrs, d, N);
T_gf = 10 * rand(nrs, N);
Tijs = 10 * rand(d, num_edges);
problem_data.Tijs = Tijs;


% Construct a manifold structure representing the product of groups of
% rotations with the Euclidean space for A. We optimize simultaneously
% for the reference cloud and for the rotations that affect each of the
% measured clouds. Notice that there is a group invariance because
% there is no way of telling which orientation the reference cloud
% should be in.
tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
M = productmanifold(tuple);

% Setup the problem structure with manifold M and cost+grad functions.
problem.M = M;
problem_data.sz = [nrs, d, N];
% problem_struct = problem
problem.cost = @(x) cost(x, problem_data);
problem.grad = @(x) grad(x, problem_data);
problem.hess = @(x, u) hess(x,u, problem_data);

% For debugging, it's always nice to check the gradient a few times.
disp("Checking gradient...");
figure(1);
checkgradient(problem);
% pause;
disp("Checking hessian...");
figure(2);
checkhessian(problem);
% pause;

% Call a solver on our problem. This can probably be much improved if a
% clever initial guess is used instead of a random one.
X = trustregions(problem);
T = X.T;
R = X.R;



%%
% Define the cost function here. Points on the manifold M are
% structures with fields X.A and X.R, containing matrices of sizes
% respectively nxm and nxnxN. 
% The store structure (the caching system) is used to keep the residue 
% matrix E in memory, as it is also used in the computation of the 
% gradient and of the Hessian. This way, we prevent redundant computations.
function [f] = cost(X, problem_data)
%     f = rsom_cost_rot_stiefel(R, problem);
    f = rsom_cost_base(X, problem_data);
end %cost

%%
% Riemannian gradient of the cost function.
function [g] = grad(X, problem_data)
    R = X.R;
    T = X.T;
    Tijs = problem_data.Tijs;
    edges = problem_data.edges;
    [problem_structs.P, problem_structs.frct] = ...
        make_step1_p_fct(T, Tijs, edges);
    [problem_structs.LR, problem_structs.PR, problem_structs.BR] = ...
        make_LR_PR_BR_noloops(R, Tijs, edges);
    g.R = rsom_rgrad_rot_stiefel(R, problem_structs);
    g.T = rsom_rgrad_transl_stiefel(T, problem_structs);
end %grad


%%
% It is not necessary to define the Hessian of the cost. We do it
% mostly to illustrate how to do it and to study the spectrum of the
% Hessian at the solution (see further down).
function [h] = hess(X, Xdot, problem_data)
    R = X.R;
    T = X.T;
    Rdot = Xdot.R;
    Tdot = Xdot.T;
    % Careful: tangent vectors on the rotation group are represented as
    % skew symmetric matrices. To obtain the corresponding vectors in
    % the ambient space, we need a little transformation. This
    % transformation is typically not needed when we compute the
    % formulas for the gradient and the Hessian directly in Riemannian
    % form instead of resorting the egrad2rgrad and ehess2rhess. These
    % latter tools are convenient for prototyping but are not always
    % the most efficient form to execute the computations.
    Tijs = problem_data.Tijs;
    edges = problem_data.edges;
    [problem_structs.P, problem_structs.frct] = ...
        make_step1_p_fct(T, Tijs, edges);
    [problem_structs.LR, problem_structs.PR, problem_structs.BR] = ...
        make_LR_PR_BR_noloops(R, Tijs, edges);

    h.R = rsom_rhess_rot_stiefel(R, Rdot, problem_structs);
    h.T = rsom_rhess_transl_stiefel(T, Tdot, problem_structs);
end %hess






