function scale_formulation_manopt

% importmanopt;


nrs = 4;
d = 3;
N = 5;



num_edges = 8;
G = graph(true(N), 'omitselfloops'); % Alternative without self-loops
p = randperm(numedges(G), num_edges);
G = graph(G.Edges(p, :));
edges = table2array(G.Edges);


Tijs_vec = 10 * rand(d, num_edges);
problem_data.Tijs = Tijs_vec;
problem_data.edges = edges;


tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
tuple.lambda = euclideanfactory(num_edges, 1);
M = productmanifold(tuple);

problem_data.rho = 0.0;

problem.M = M;
problem_data.sz = [nrs, d, N];
problem.cost = @(x) ssom_cost(x, problem_data);
problem.egrad = @(x) ssom_egrad(x, problem_data);
% problem.grad = @(x) ssom_rgrad(x, problem_data);
% problem.hess = @(x, u) hess_genproc(x, u, problem_data);

fprintf("\n");
disp("Gradient check easy");
X_gradcheck.R = eye3d(nrs, d, N);
X_gradcheck.T = zeros(nrs, N);
X_gradcheck.lambda = ones(num_edges, 1);
checkgradient(problem, X_gradcheck);

fprintf("\n");
disp("Gradient check rand");
checkgradient(problem);
% checkhessian(problem);

% problem_data.R_gt = eye3d(nrs, d, N);
% problem_data.T_gt = zeros(nrs, N);
% problem_data.lambda_gt = zeros(num_edges, 1);



% trustregions(problem);


end %file function



