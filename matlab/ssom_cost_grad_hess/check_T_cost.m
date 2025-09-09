function check_T_cost

nrs = 10;
d = 3;
N = 7;

mindeg = 3;

testdata = testNetwork_params(3, N, 'banded', mindeg); %4 would be the default



num_edges = size(testdata.E, 1);
% tijs = 5 * rand(num_edges, 3);

R = make_rand_stiefel_3d_array(nrs, d, N);
T = 20 * rand(nrs, N);
lambdas = 10 * rand(num_edges, 1);

X.R = R;
X.T = T;
X.lambda = lambdas;

%% setup problem_data from testdata (i.e., testNetwork output)
testdata.rho = 0.0;
problem_data = testdata;
problem_data.edges = testdata.E;
problem_data.tijs = G2T(testdata.gij);
%%
lhs = ssom_cost(X, problem_data);

[L, P, B] = ssom_cost_T_reformulation(R, lambdas, problem_data);

rhs = trace(T * L * T') + 2 * trace(T * P') + trace(B); % !! P'

disp("[lhs, rhs]")
disp([lhs, rhs])



end %file function


function [L, P, B] = ssom_cost_T_reformulation(R, lambdas, problem_data)


% g = zeros(size(T));

N = size(R, 3);
nrs = size(R, 1);

L = zeros(N,N);
P = zeros(nrs,N);
% BR_const = zeros(d,d);
B = 0;


num_edges = size(problem_data.edges,1);
for e = 1:num_edges
    ii = problem_data.edges(e,1);
    jj = problem_data.edges(e,2);
    bij = zeros(N,1);
    bij(ii, 1) = 1;
    bij(jj, 1) = -1;
    tij = problem_data.tijs(:, e);
    lambda_e = lambdas(e, 1);
    L = L + (bij * bij');
    Ri = R(:,:,ii);
    P = P + lambda_e * (Ri * tij * bij');
    % B = B + lambda_e ^ 2 * (Ri * tij * tij' * Ri');
    B = B + lambda_e ^ 2 * (tij * tij');
end


end