function test_check_ssom_test_ehess_T_lambda_v1()

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

gradf_T= @(t) egrad_T(R0, T0, lambda(t), problem_data);
hessf_T_lambda = @(t) ssom_ehess_t_lambda(R0, T0, lambda(t), dLambda(t), problem_data);

funCheckDer(gradf_T, hessf_T_lambda)

end % file function

function g=egrad_T(R, T, lambdas, problem_data)
% g=x*(problem_data.LR+problem_data.LR')+(problem_data.PR)';
tijs_scaled = make_tijs_scaled(lambdas, problem_data.tijs);
[LR, PR] = make_LR_PR_BR_noloops(R, tijs_scaled, problem_data.edges);
g=T*(LR+LR')+PR';
end

function h = ssom_ehess_t_lambda(R, T, ~, lambdas_dot, problem_data)
edges = problem_data.edges;

h = zeros(size(T));
N = size(R,3);
num_edges = size(edges, 1);
for e = 1:num_edges
    ii = edges(e,1);
    jj = edges(e,2);
    Ri = R(:,:,ii);
    BIJ = zeros(N,1);
    BIJ(ii) = 1;
    BIJ(jj) = -1;
    %
    tij = problem_data.tijs(:, e);
    w_ij = BIJ * lambdas_dot(e) * tij' * Ri';
    h = h +  2 * w_ij'; % w_ij transpose!
end
end

