function scale_formulation_relu

d = 3;
nrs = 4;
N = 5;

%graph random init
num_edges = 8;
G = graph(true(N), 'omitselfloops'); % Alternative without self-loops
% resetRands(0);
p = randperm(numedges(G), num_edges);
G = graph(G.Edges(p, :));
edges = table2array(G.Edges);


% resetRands(0);
lambdas = 5 * rand(num_edges, 1);
% resetRands(0);
Tijs_vec = 5 * ones(d, num_edges) - 10 * rand(d, num_edges);

% resetRands(0);
R_globalframe = make_rand_stiefel_3d_array(nrs, d, N);
% resetRands(0);
T_globalframe = 10 * rand(nrs, N);

% resetRands(0);
rho = 1000 * rand(1,1);

% setup X
X.R = R_globalframe;
X.T = T_globalframe;
X.lambda = lambdas;

% setup problem_data
problem_data.edges = edges;
problem_data.Tijs = Tijs_vec;
problem_data.rho = rho;


cost_lambda = ssom_cost_norms(X, problem_data);

X.lambda = -lambdas;
cost_minus_lambda = ssom_cost_norms(X, problem_data);

X.lambda = zeros(size(X.lambda));

cost_lambda_0 = ssom_cost_norms(X, problem_data);

%compute cost by summing norms (prior to trace reformulation)
% cost_norm = ssom_cost_norms(X, problem_data);

% disp("cost_norm")
% disp(cost_norm)

disp("cost_lambda")
disp(cost_lambda)

disp("cost_minus_lambda")
disp(cost_minus_lambda)

disp("cost_lambda > cost_minus_lambda")
disp(cost_lambda > cost_minus_lambda)

disp("cost_lambda_0")
disp(cost_lambda_0)

disp("cost_lambda > cost_lambda_0")
disp(cost_lambda > cost_lambda_0)

end %file function
