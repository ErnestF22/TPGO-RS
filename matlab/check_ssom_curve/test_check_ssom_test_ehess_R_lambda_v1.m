function test_check_ssom_test_ehess_R_lambda_v1()

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

problem_data = struct('sz', sz, 'edges', edges, 'tijs', tijs, 'rho', rho);

e = size(problem_data.edges,1);

%%

lambda0=rand(e,1,1);
vLambda0=rand(e,1,1);

R0 = make_rand_stiefel_3d_array(nrs,d,N);

T0 = rand(nrs,N);

[lambda,dLambda,~,~,~]=real_geodFun(lambda0, vLambda0);

gradf_R = @(t) matStackH(egrad_R(R0, T0, lambda(t), problem_data));
hessf_R_lambda = @(t) matStackH(ssom_ehess_R_lambda(R0, T0, lambda(t), dLambda(t), problem_data));

funCheckDer(gradf_R, hessf_R_lambda)

end

function g = egrad_R(~, T, lambdas, problem_data)
nrs = size(T, 1);
d = size(problem_data.tijs, 1);
N = size(T, 2);

P = zeros(nrs, d*N);

tijs_scaled = make_tijs_scaled(lambdas, problem_data.tijs);

idx_col_p = reshape(1:d*N, [], N)';

num_edges = size(problem_data.edges,1);
for e = 1:num_edges
    ii = problem_data.edges(e,1);
    jj = problem_data.edges(e,2); 
    T_j = T(:, jj);
    T_i = T(:, ii);
    tij = tijs_scaled(:,e);
    P_e = 2 * (T_i * tij' - T_j * tij');
    P(:, idx_col_p(ii, :)) = ...
        P(:, idx_col_p(ii, :)) + P_e;
end

g=matUnstackH(P,d);
end

function h = ssom_ehess_R_lambda(~, T, ~, lambdas_dot, problem_data)
nrs = size(T, 1);
d = size(problem_data.tijs, 1);
N = size(T, 2);

Ph = zeros(nrs, d*N);

tijs_dot_scaled = make_tijs_scaled(lambdas_dot, problem_data.tijs);

idx_col_p = reshape(1:d*N, [], N)';

num_edges = size(problem_data.edges,1);
for e = 1:num_edges
    ii = problem_data.edges(e,1);
    jj = problem_data.edges(e,2); 
    T_j = T(:, jj);
    T_i = T(:, ii);
    tij_dot = tijs_dot_scaled(:,e);
    P_e = 2 * (T_i * tij_dot' - T_j * tij_dot');
    Ph(:, idx_col_p(ii, :)) = ...
        Ph(:, idx_col_p(ii, :)) + P_e;
end

h=matUnstackH(Ph,d);
end