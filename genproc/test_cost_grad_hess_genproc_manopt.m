close all
clear
clc

resetRands(2);

d = 3;
nrs = 8;
N = 80;

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
problem.cost = @(x) cost_genproc(x, problem_data);
problem.grad = @(x) grad_genproc(x, problem_data);
problem.hess = @(x, u) hess_genproc(x, u, problem_data);

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
% X = trustregions(problem);
% T = X.T;
% R = X.R;













