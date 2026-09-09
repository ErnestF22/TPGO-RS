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

ssom_cost_qp2 = ssom_cost_qp(X, problem_data);
disp("ssom_cost_qp2")
disp(ssom_cost_qp2)

ssom_cost_full2 = ssom_cost_full(X, problem_data);
disp("ssom_cost_full2")
disp(ssom_cost_full2)

end %file function

function [A, B] = make_a_b_qp(R, problem_data, ii, jj, ee)

num_edges = size(problem_data.edges, 1);
N = size(R, 3);
A = zeros(3, 3*N);
B = zeros(3, num_edges);

tijs = problem_data.tijs;

A(:, (ii-1)*3+1:ii*3) = eye(3);
A(:, (jj-1)*3+1:jj*3) = -eye(3);

B(:, ee) = R(:,:,ii) * tijs(:,ee);

end

function cost_out = ssom_cost_qp(X, problem_data)
num_edges = size(problem_data.edges, 1);
edges = problem_data.edges;
R = X.R;
T = X.T;
lambdas = X.lambda;
cost_out = 0.0;

for ee = 1:num_edges
    ii = edges(ee, 1);
    jj = edges(ee, 2);

    [A_ee, B_ee] = make_a_b_qp(R, problem_data, ii, jj, ee);

    cost_increment = norm([A_ee, B_ee] * [T(:); lambdas(:)]);
    cost_out = cost_out + cost_increment*cost_increment;

end

end


function cost_out = ssom_cost_full(X, problem_data)
num_edges = size(problem_data.edges, 1);
edges = problem_data.edges;
R = X.R;
T = X.T;
lambdas = X.lambda;
% cost_out = 0.0;

A_ee_full = [];
B_ee_full = [];

for ee = 1:num_edges
    ii = edges(ee, 1);
    jj = edges(ee, 2);

    [A_ee, B_ee] = make_a_b_qp(R, problem_data, ii, jj, ee);

    A_ee_full = [A_ee_full; A_ee];
    B_ee_full = [B_ee_full; B_ee];

end


cost_out = norm([A_ee_full, B_ee_full] * [T(:); lambdas(:)]); % Accumulate the cost

cost_out = cost_out * cost_out;

end