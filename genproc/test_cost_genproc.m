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
problem.edges = edges;
% problem.num_edges = num_edges;

R_gf = make_rand_stiefel_3d_array(nrs, d, N);
T_gf = 10 * rand(nrs, N);
Tijs = 10 * rand(d, num_edges);
% problem.R_gf = R_gf;
% problem.T_gf = T_gf;
problem.Tijs = Tijs;


% Construct a manifold structure representing the product of groups of
% rotations with the Euclidean space for A. We optimize simultaneously
% for the reference cloud and for the rotations that affect each of the
% measured clouds. Notice that there is a group invariance because
% there is no way of telling which orientation the reference cloud
% should be in.
tuple.R = stiefelfactory(nrs, d, N);
tuple.A = euclideanfactory(nrs, N);
M = productmanifold(tuple);

% Setup the problem structure with manifold M and cost+grad functions.
problem.M = M;
problem.sz = [nrs, d, N];
% problem_struct = problem
[problem.P, problem.frct] = make_step1_p_fct(T_gf, Tijs, edges);
[problem.LR, problem.PR, problem.BR] = make_LR_PR_BR_noloops(R_gf, Tijs, edges);
problem.cost = @(x, store) cost(x,problem,store);
problem.grad = @(x, store) grad(x,problem,store);
problem.hess = @(x, u, store) hess(x,u,problem,store);

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
A = X.A;
R = X.R;



%%
% Define the cost function here. Points on the manifold M are
% structures with fields X.A and X.R, containing matrices of sizes
% respectively nxm and nxnxN. The store structure (the caching system)
% is used to keep the residue matrix E in memory, as it is also used in
% the computation of the gradient and of the Hessian. This way, we
% prevent redundant computations.
function [f, store] = cost(X, problem, store)
    R = X.R;
    A = X.A;
    Tijs = problem.Tijs;
    edges = problem.edges;
%     f = rsom_cost_rot_stiefel(R, problem);
    f = rsom_cost_base(R, A, Tijs, edges);
end %cost

%%
% Riemannian gradient of the cost function.
function [g, store] = grad(X, problem, store)
    R = X.R;
    A = X.A;
    g.R = rsom_rgrad_rot_stiefel(R, problem);
    g.A = rsom_rgrad_transl_stiefel(A, problem);
end %grad


%%
% It is not necessary to define the Hessian of the cost. We do it
% mostly to illustrate how to do it and to study the spectrum of the
% Hessian at the solution (see further down).
function [h, store] = hess(X, Xdot, problem, store)
    R = X.R;
    A = X.A;
    Rdot = Xdot.R;
    Adot = Xdot.A;
    % Careful: tangent vectors on the rotation group are represented as
    % skew symmetric matrices. To obtain the corresponding vectors in
    % the ambient space, we need a little transformation. This
    % transformation is typically not needed when we compute the
    % formulas for the gradient and the Hessian directly in Riemannian
    % form instead of resorting the egrad2rgrad and ehess2rhess. These
    % latter tools are convenient for prototyping but are not always
    % the most efficient form to execute the computations.
    h.R = rsom_rhess_rot_stiefel(R, Rdot, problem);
    h.A = rsom_rhess_transl_stiefel(A, Adot, problem);
end %hess






