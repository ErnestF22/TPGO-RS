function test_qp_reformulation_T_lambda

sigma_noise = 0.0;

N = 5;
mindeg = 2;
testdata = testNetwork_params(3, N, 'banded', mindeg, 0.0, sigma_noise);
testdata.mindeg = mindeg;
testdata.sigma = sigma_noise;
testdata.rho = 0; % compensation part
testdata.a = 2.0; % log scale-compensation param

X.R = G2R(testdata.gi);
X.T = G2T(testdata.gi);
X.lambda = testdata.lambdaij;

resetRands(posixtime(datetime("now")));

X.R = rand(size(X.R));
X.T = rand(size(X.T));
X.lambda = rand(size(X.lambda));

problem_data.tijs = G2T(testdata.gij);
problem_data.tijs = rand(size(problem_data.tijs));
problem_data.edges = testdata.E;
problem_data.rho = 0.0;
problem_data.a = testdata.a;

ssom_cost_orig = ssom_cost(X, problem_data);
disp("ssom_cost_orig")
disp(ssom_cost_orig)

ssom_cost_qp2 = ssom_T_lambda_cost_qp_sum(X, problem_data);
disp("ssom_cost_qp2")
disp(ssom_cost_qp2)

ssom_cost_full2 = ssom_T_lambda_cost_qp_matricial(X, problem_data);
disp("ssom_cost_full2")
disp(ssom_cost_full2)

end %file function

