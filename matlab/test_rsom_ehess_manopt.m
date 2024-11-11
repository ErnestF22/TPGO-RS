clc;
clear;
close all;

nrs = 4;
d = 3;
N = 5;



num_edges = 8;
G = graph(true(N), 'omitselfloops'); % Alternative without self-loops
p = randperm(numedges(G), num_edges);
G = graph(G.Edges(p, :));
edges = table2array(G.Edges);


Tijs = 10 * rand(d, num_edges);
problem_data.Tijs = Tijs;
problem_data.edges = edges;




problem_data.Tijs = Tijs;
problem_data.edges = edges;




tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
M = productmanifold(tuple);

% Setup the problem structure with manifold M and cost+grad functions.
problem.M = M;
problem_data.sz = [nrs, d, N];
problem.cost = @(x) cost_genproc(x, problem_data);
problem.egrad = @(x) egrad_genproc(x, problem_data);
problem.ehess = @(x, u) ehess_genproc(x, u, problem_data);



figure(1)
checkgradient(problem);

figure(2)
checkhessian(problem);


% X = trustregions(problem);
% T_manopt_out = X.T;
% R_manopt_out = X.R;




