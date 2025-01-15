function test_check_ssom_test_ehess_lambda_T_v1

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

N = problem_data.sz(3);
% nrs = problem.sz(1);
d = problem_data.sz(2);
e = size(problem_data.edges,1);

lambda0=rand(e,1,1);
% vLambda0=rand(e,1,1);

R0 = eye3d(d,d,N); %first d with nrs later
% vR0 = eye3d(d,d,N);

T0 = rand(d,N);
vT0 = rand(d,N);

[T,dT,~,~,~]=real_geodFun(T0, vT0);

gradf_T= @(t) egrad_lambda(R0, T(t), lambda0, problem_data);
hessf_T_lambda = @(t) ssom_ehess_lambda_t(R0, T(t), dT(t), lambda0, problem_data);

funCheckDer(gradf_T, hessf_T_lambda)

end %file function

function g_lambda = egrad_lambda(R, T, lambdas, problem_data)
    edges = problem_data.edges;
    tijs = problem_data.tijs;
    rho = problem_data.rho;


    g_lambda = zeros(length(lambdas), 1);
    
    num_edges = size(edges, 1);
    for ee = 1:num_edges
        ii = edges(ee, 1);
        jj = edges(ee, 2);
        lambda_e = lambdas(ee);
        tij_e = tijs(:, ee);
        T_i = T(:, ii);
        T_j = T(:, jj);
        R_i = R(:, :, ii);
        a = T_i - T_j;
        b = R_i * tij_e;
        base_part = 2*(b' * b * lambda_e + a' * b);
        relu_part = 0.0;
        if (ssom_relu_argument(lambda_e)>0)
            relu_part = -1.0;
        end
        g_lambda(ee) = base_part + rho * relu_part;
    end
end

function h = ssom_ehess_lambda_t(R, ~, dT, lambdas, problem_data)
edges = problem_data.edges;

h = zeros(size(lambdas));
num_edges = size(edges, 1);
    for ee = 1:num_edges
        ii = edges(ee, 1);
        jj = edges(ee, 2);
        % lambda_e = lambdas(ee);
        tij_e = problem_data.tijs(:, ee);
        dT_i = dT(:, ii);
        dT_j = dT(:, jj);
        R_i = R(:, :, ii);
        adot = dT_i - dT_j;
        b = R_i * tij_e;
        h_ee = 2*(adot' * b);
        h(ee) = h_ee;
    end
end