function ssom_fix_cpp_grad_hess

startXvec = - 2 * load("../matlab/data/ssom_testdata_noisy/harder/tdata_n5_mindeg2_sigma00/ssom_x_start.csv");

nrs = 3;
d = 3;
N = 5;

mindeg = 2;

rho = 5.0; %OBS. No rho to be loaded!


problem_data = testNetwork_params(3, N, 'banded', mindeg);
problem_data.edges = problem_data.E;
problem_data.tijs = G2T(problem_data.gij);
problem_data.rho = rho;

num_edges = size(problem_data.E, 1);

tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
tuple.lambda = euclideanfactory(num_edges, 1);
M = productmanifold(tuple);


% Setup the problem structure with manifold M and cost+grad functions.
problem.M = M;
problem.cost = @(x) ssom_cost(x, problem_data);
% problem.egrad = @(x) ssom_egrad(x, problem_data);
problem.grad = @(x) ssom_rgrad(x, problem_data);
% problem.ehess = @(x, u) ssom_ehess_genproc(x, u, problem_data);
problem.hess = @(x, u) ssom_rhess_genproc(x, u, problem_data);


startX = convertXtoRTLambdas(startXvec, nrs, d, N);

% randU = problem.M.randvec(startX);
% randUvec = vectorizeXrtlambdas(randU);
% save("../matlab/data/ssom_testdata_noisy/harder/tdata_n5_mindeg2_sigma00/randUvec.csv", "randUvec", "-ascii");
load("../matlab/data/ssom_testdata_noisy/harder/tdata_n5_mindeg2_sigma00/randUvec.csv", "randUvec");

randU = convertXtoRTLambdas(randUvec, nrs, d, N);

disp("problem.cost(startX)");
disp(problem.cost(startX));

g = problem.grad(startX);
disp("g.R");
disp(g.R);
disp("g.T");
disp(g.T);
disp("g.lambda");
disp(g.lambda);

h = problem.hess(startX, randU);
disp("h.R");
disp(h.R);
disp("h.T");
disp(h.T);
disp("h.lambda");
disp(h.lambda);


problem_data_next = problem_data;
problem_data_next.sz = [d + 1, d, N];
ssom_pim_hessian_genproc(startX, problem_data_next, 1e-5, 5000)

end % file function
