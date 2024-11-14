function scale_formulation

% resetRands(2);

d = 3;
nrs = 4;
N = 5;

%graph random init
num_edges = 8;
G = graph(true(N), 'omitselfloops'); % Alternative without self-loops
p = randperm(numedges(G), num_edges);
G = graph(G.Edges(p, :));
edges = table2array(G.Edges);


lambdas = 10 * rand(num_edges, 1);
Tijs_vec = 10 * rand(d, num_edges);

R_globalframe = make_rand_stiefel_3d_array(nrs, d, N);
T_globalframe = 10 * rand(nrs, N);

rho = 100 * rand(1,1);% setup X

% setup X
X.R = R_globalframe;
X.T = T_globalframe;
X.lambda = lambdas;

% setup problem_data
problem_data.edges = edges;
problem_data.Tijs = Tijs_vec;
problem_data.rho = rho;


X.R = R_globalframe;
X.T = T_globalframe;
X.lambda = lambdas;

% setup problem_data
problem_data.edges = edges;
problem_data.Tijs = Tijs_vec;
problem_data.rho = rho;

cost_norms = ssom_cost_norms(X, problem_data);

cost_trace = ssom_cost(X, problem_data);

cost_lambda = ssom_cost_lambda(X, problem_data);

disp("cost_norms")
disp(cost_norms)

disp("cost_trace")
disp(cost_trace)

disp("cost_lambda")
disp(cost_lambda)


