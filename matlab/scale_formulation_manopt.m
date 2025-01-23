function scale_formulation_manopt

% importmanopt;


nrs = 4;
d = 3;
N = 5;

sz=[nrs,d,N];

%graph random init
num_edges = 8;
G = graph(true(N), 'omitselfloops'); % Alternative without self-loops
p = randperm(numedges(G), num_edges);
G = graph(G.Edges(p, :));
edges = table2array(G.Edges);

tijs = 10 * rand(d, num_edges);


rho = 1.0; %TODO: make this rand() later

% variables random generation/init
tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
tuple.lambda = euclideanfactory(num_edges, 1);
M = productmanifold(tuple);
% problem.M = M;
% X = M.rand();

problem_data = struct('sz', sz, 'edges', edges, 'tijs', tijs);

%g.lambda
% [aL, bL, cL] = makeABClambda(X, problem_data);
% problem_data.aL = aL; problem_data.bL = bL; problem_data.cL = cL;

problem_data.rho = rho; %ReLU() part should not be needed for Hessian tests

%% output (problem_curve_data) definition

problem.M = M;
problem_data.sz = [nrs, d, N];
problem.cost = @(x) ssom_cost(x, problem_data);
% problem.egrad = @(x) ssom_egrad(x, problem_data);
% problem.ehess = @(x, u) ssom_ehess_genproc(x, u, problem_data);
problem.grad = @(x) ssom_rgrad(x, problem_data);
problem.hess = @(x, u) ssom_rhess_genproc(x, u, problem_data);

fprintf("\n");
disp("Gradient check easy");
X_gradcheck.R = eye3d(nrs, d, N);
X_gradcheck.T = ones(nrs, N);
X_gradcheck.lambda = zeros(num_edges, 1);
figure(1)
checkgradient(problem, X_gradcheck);

fprintf("\n");
disp("Gradient check rand");
figure(2)
checkgradient(problem);

% close all;

fprintf("\n");
disp("Hessian check easy");
% X_gradcheck.R = eye3d(nrs, d, N);
% X_gradcheck.T = zeros(nrs, N);
% X_gradcheck.lambda = zeros(num_edges, 1);
% Xdot_gradcheck.R = eye3d(nrs, d, N);
% Xdot_gradcheck.T = zeros(nrs, N);
% Xdot_gradcheck.lambda = zeros(num_edges, 1);

figure(3)
X_hesscheck = M.rand();
Xdot_hesscheck = M.rand();
X_hesscheck.R = eye3d(nrs, d, N);
Xdot_hesscheck.R = stiefel_randTangentNormVector(X_hesscheck.R);
% only case this works is with the following T, lambda and R = zeros()
X_hesscheck.T = ones(nrs, N);
X_hesscheck.lambda = 2 * ones(num_edges, 1);
checkhessian(problem, X_hesscheck, Xdot_hesscheck);


fprintf("\n");
figure(4)
disp("Hessian check rand");
checkhessian(problem);

% problem_data.R_gt = eye3d(nrs, d, N);
% problem_data.T_gt = zeros(nrs, N);
% problem_data.lambda_gt = zeros(num_edges, 1);



% trustregions(problem);


end %file function

