function test_check_ssom_test_ehess_lambda_R_v1

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

rho = 2.0; %TODO: make this rand() later

problem_data = struct('sz', sz, 'edges', edges, 'tijs', tijs, 'rho', rho);

e = size(problem_data.edges,1);

%%

lambda0=rand(e,1,1);
% vLambda0=rand(e,1,1);

R0 = randrot(d,N);
vR0 = rand(d,d,N);

T0 = rand(d,N);
% vT0 = rand(d,N);

[R,dR,~,~,~]=real_geodFun(R0, vR0);

gradf_lambda= @(t) egrad_lambda(R(t), T0, lambda0, problem_data);
hessf_lambda_R = @(t) ssom_ehess_lambda_R(R(t), dR(t), T0, lambda0, problem_data);

funCheckDer(gradf_lambda, hessf_lambda_R)

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

function h = ssom_ehess_lambda_R(R, dR, T, lambdas, problem_data)
edges = problem_data.edges;

h = zeros(size(lambdas));
num_edges = size(edges, 1);
    for ee = 1:num_edges
        ii = edges(ee, 1);
        jj = edges(ee, 2);
        lambda_e = lambdas(ee);
        % lambda_e_sq = lambda_e * lambda_e;
        tij_e = problem_data.tijs(:, ee);
        T_i = T(:, ii);
        T_j = T(:, jj);
        dR_i = dR(:, :, ii);
        R_i = R(:, :, ii);
        a = T_i - T_j;
        b = R_i * tij_e;
        bdot = dR_i * tij_e;
        h_ee_half = lambda_e * (bdot' * b) + lambda_e * (b' * bdot)  + ...
                    a' * bdot;
        h(ee) = 2 * h_ee_half;
    end
    % 2λij ˙bT b + 2λij bT ˙b + 2λij aT ˙b
end