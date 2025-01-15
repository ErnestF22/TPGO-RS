function test_check_ssom_test_ehess_T_lambda_v1()
problem=test_check_ssom();
% curve=test_check_ssom_curve(problem);

% clean problem struct of unwanted members
problem_data.sz = problem.sz;
problem_data.edges = problem.edges;
problem_data.tijs = problem.tijs;

problem = [];

N = problem_data.sz(3);
% nrs = problem.sz(1);
d = problem_data.sz(2);
e = size(problem_data.edges,1);

lambda0=rand(e,1,1);
vLambda0=rand(e,1,1);

R0 = eye3d(d,d,N); %first d with nrs later
% vR0 = eye3d(d,d,N);

T0 = rand(d,N);
% vT0 = rand(d,N);

[lambda,dLambda,~,~,~]=real_geodFun(lambda0, vLambda0);

% curve.c=@(t) T(t);
% curve.dc=@(t) dT(t);
% curve.ddc=@(t) ddT(t);

% % f=@(t) problem.cost(curve.c(t));
% gradf=@(t) egrad_T(T0);
% df=@(t) sum(stiefel_metric([],gradf(t),curve.dc(t)));
% % funCheckDer(f,df)
% ehessf = @(t) ssom_ehess_t_lambda(R0, T0, dLambda(t), ddLambda(t));
% ddf_1 = @(t) stiefel_metric([], ehessf(t), curve.dc(t), 'euclidean');
% ddf_2 = @(t) stiefel_metric([], gradf(t), curve.ddc(t), 'euclidean');
% ddf = @(t) sum(ddf_1(t) + ddf_2(t));    
% 
% funCheckDer(df,ddf)

%% version with fixed T, geodesic lambda

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

