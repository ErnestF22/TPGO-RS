function test_solve_ssom_T_lambdas_qp

sigma_noise = 0.0;

N = 5;
mindeg = 2;
testdata = testNetwork_params(3, N, 'banded', mindeg, 0.0, sigma_noise);
testdata.mindeg = mindeg;
testdata.sigma = sigma_noise;
testdata.rho = 0; % compensation part
testdata.a = 2.0; % log scale-compensation param

num_edges = size(testdata.E, 1);


R = G2R(testdata.gi);
% X.T = G2T(testdata.gi);
% X.lambda = testdata.lambdaij;

% resetRands(posixtime(datetime("now")));
% 
% X.R = rand(size(X.R));
% X.T = rand(size(X.T));
% X.lambda = rand(size(X.lambda));


problem_data.tijs = G2T(testdata.gij);
problem_data.tijs = rand(size(problem_data.tijs));
problem_data.edges = testdata.E;
problem_data.rho = 0.0;
problem_data.a = testdata.a;

problem_data.N = N;
problem_data.num_edges = num_edges;
[T, lambdas] = solve_ssom_T_lambdas_qp(R, problem_data);

X_out.R = R; % unchanged
X_out.T = T;
X_out.lambda = lambdas;

disp("ssom_cost(X_out, problem_data)")
disp(ssom_cost(X_out, problem_data))

% resetRands(posixtime(datetime("now")));
% X_out.T = rand(size(T));
% 
% disp("ssom_cost(X_out, problem_data)")
% disp(ssom_cost(X_out, problem_data))

end %file function





