function test_check_ssom_test_ehess_T_R_v1

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
% vLambda0=rand(e,1,1);

R0 = make_rand_stiefel_3d_array(nrs,d,N);
vR0 = rand(nrs, d, N);

T0 = rand(nrs,N);

[R,dR,~,~,~]=real_geodFun(R0, vR0);

gradf_T= @(t) egrad_T(R(t), T0, lambda0, problem_data);
hessf_T_R = @(t) ssom_ehess_T_R(R(t), dR(t), T0, lambda0, problem_data);

funCheckDer(gradf_T, hessf_T_R)

end % file function

function g=egrad_T(R, T, lambdas, problem_data)
% g=x*(problem_data.LR+problem_data.LR')+(problem_data.PR)';
tijs_scaled = make_tijs_scaled(lambdas, problem_data.tijs);
[LR, PR] = make_LR_PR_BR_noloops(R, tijs_scaled, problem_data.edges);
g=T*(LR+LR')+PR';
end

function h = ssom_ehess_T_R(~, dR, ~, lambdas, problem_data)
tijs_scaled = make_tijs_scaled(lambdas, problem_data.tijs);
[~, PR_dot] = make_LR_PR_BR_noloops(dR, tijs_scaled, problem_data.edges);
h=PR_dot';
end

