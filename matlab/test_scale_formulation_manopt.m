function test_scale_formulation_manopt

% importmanopt;


%% testNetwork

N = 5;
mindeg = 3;
testdata = testNetwork_params(3, N, 'banded', mindeg);
testdata.mindeg = 3;

X_gt.lambda = testdata.lambdaijtruth;
X_gt.R = G2R(testdata.gitruth);
X_gt.T = G2T(testdata.gitruth);

disp("X_gt")
disp(X_gt)

tijs = G2T(testdata.gij);
edges = testdata.E;


%%

nrs = size(X_gt.T, 1);
d = size(tijs, 1); % d = 3
N = size(X_gt.R, 3);

sz=[nrs,d,N];

num_edges = size(edges, 1);


rho = 1.0; %TODO: set this properly

% variables random generation/init
tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
tuple.lambda = euclideanfactory(num_edges, 1);
M = productmanifold(tuple);

problem_data = struct('sz', sz, 'edges', edges, 'tijs', tijs);

problem_data.rho = rho; 

%% check gradient


problem.M = M;
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

close all;

%% check hessian

fprintf("\n");
disp("Hessian check easy");
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
disp("Hessian check rand");
checkhessian(problem);

close all;

problem_data.R_gt = G2R(testdata.gitruth);
problem_data.T_gt = G2T(testdata.gitruth);
problem_data.lambda_gt = testdata.lambdaijtruth;

disp("problem_data.lambda_gt")
disp(problem_data.lambda_gt)


%% call RS-based optimization function

ssom_genproc(problem_data)




end %file function
