function test_check_ssom_test_ehess_T_lambda_v1()


nrs = 3;
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

lambdas = 10 * rand(num_edges, 1);

rho = 1.0; %TODO: make this rand() later

tijs_scaled = make_tijs_scaled(lambdas, tijs); %!!
problem_data = struct('sz', sz, 'edges', edges, 'tijs', tijs, 'Tijs', tijs_scaled);

% variables random generation/init
tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
tuple.lambda = euclideanfactory(num_edges, 1);
M = productmanifold(tuple);
% problem.M = M;
x = M.rand();
lambdas = x.lambda;
T = x.T;
R = x.R;


%g.R
[problem_data.P, problem_data.frct] = ...
    make_step1_p_fct(T, tijs_scaled, edges);
%g.T
problem_data.Tijs = tijs_scaled;
[problem_data.LR, problem_data.PR, problem_data.BR] = ...
    make_LR_PR_BR_noloops(R, tijs_scaled, edges);
%g.lambda
problem_data.T = T;
problem_data.R = R;
[aL, bL, cL] = makeABClambda(x, problem_data);
problem_data.aL = aL; problem_data.bL = bL; problem_data.cL = cL;


gT = rsom_egrad_transl_stiefel(T, problem_data);


vLambda0 = rand(size(lambdas));
[lambda,dLambda,~,~,~]=real_geodFun(lambdas, vLambda0);


LR = zeros(N,N);
PR = zeros(N,nrs);
BR_const = zeros(d,d);
for ee = 1:num_edges
    BIJ = zeros(N,1);
    ii = edges(ee, 1);
    jj = edges(ee, 2);
    BIJ(ii) = 1;
    BIJ(jj) = -1;
    Tij = problem_data.Tijs(:, ee);
    Ri = R(:,:,ii);
    %LR
    LR = LR + BIJ * BIJ';
    %PR
    PR = PR + 2 * BIJ * Tij' * Ri';
    %BR
    BR_const = BR_const + Tij * Tij';
end

scalar_i = 1;
scalar_j = 1;

val_cmp = gT(scalar_i, scalar_j);
val = 0;
for ee = 1:num_edges
    % ii = edges(ee, 1);
    % jj = edges(ee, 2);
    BIJ = zeros(N,1);
    ii = edges(ee, 1);
    jj = edges(ee, 2);
    BIJ(ii) = 1;
    BIJ(jj) = -1;
    %%%
    bij = BIJ(scalar_i, :); % bij is a scalar
    lambdaij = lambdas; % !! lambdas is transposed in this context
    tij = problem_data.tijs(scalar_j, :);

    val = val + bij * lambdaij' * tij';
end

disp("[val_cmp, val]")
disp([val_cmp, val])

%%
problem=test_check_ssom();
% curve=test_check_ssom_curve(problem);

N = problem.sz(3);
% nrs = problem.sz(1);
d = problem.sz(2);
e = size(problem.edges,1);

lambda0=rand(e,1,1);
vLambda0=rand(e,1,1);

R0 = eye3d(d,d,N); %first d with nrs later
vR0 = eye3d(d,d,N);

T0 = rand(d,N);
vT0 = rand(d,N);

[lambda,dLambda,~,~,ddLambda]=real_geodFun(lambda0, vLambda0);
[T,dT,~,~,ddT]=real_geodFun(T0, vT0);
[R,~,~,~,~]=rot_geodFun(R0, vR0);

curve.c=@(t) T(t);
curve.dc=@(t) dT(t);
curve.ddc=@(t) ddT(t);

% f=@(t) problem.cost(curve.c(t));
gradf=@(t) problem.egrad_T(curve.c(t));
df=@(t) sum(stiefel_metric([],gradf(t),curve.dc(t)));
% funCheckDer(f,df)
ehessf = @(t) problem.ssom_ehess_t_lambda(curve.c(t),curve.dc(t),R0, ddLambda(t));
ddf_1 = @(t) stiefel_metric([], ehessf(t), curve.dc(t), 'euclidean');
ddf_2 = @(t) stiefel_metric([], gradf(t), curve.ddc(t), 'euclidean');
ddf = @(t) sum(ddf_1(t) + ddf_2(t));    

funCheckDer(df,ddf)

%% version with fixed R0, vLambda0
% problem=test_check_ssom();
% % curve=test_check_ssom_curve(problem);
% 
% N = problem.sz(3);
% % nrs = problem.sz(1);
% d = problem.sz(2);
% e = size(problem.edges,1);
% 
% lambda0=rand(e,1,1);
% vLambda0=rand(e,1,1);
% vLambda0 = vLambda0 / norm(vLambda0);
% 
% R0 = eye3d(d,d,N); %first d with nrs later
% vR0 = rot_randTangentNormVector(R0);
% 
% T0 = rand(d,N);
% vT0 = rand(d,N);
% vT0 = vT0 / norm(vT0);
% 
% [~,~,~,~,~]=real_geodFun(lambda0, vLambda0);
% [T,dT,~,~,ddT]=real_geodFun(T0, vT0);
% [~,~,~,~,~]=rot_geodFun(R0, vR0);
% 
% curve.c=@(t) T(t);
% curve.dc=@(t) dT(t);
% curve.ddc=@(t) ddT(t);
% 
% % f=@(t) problem.cost(curve.c(t));
% egradf=@(t) problem.egrad_T(curve.c(t),R0,lambda0);
% df=@(t) sum(stiefel_metric([],egradf(t),curve.dc(t)));
% % funCheckDer(f,df)
% ehessf = @(t) problem.ssom_ehess_t_lambda(curve.c(t),curve.dc(t),R0, vLambda0);
% ddf_1 = @(t) stiefel_metric([], ehessf(t), curve.dc(t), 'euclidean');
% ddf_2 = @(t) stiefel_metric([], egradf(t), curve.ddc(t), 'euclidean');
% ddf = @(t) sum(ddf_1(t) + ddf_2(t));    
% 
% funCheckDer(df,ddf)



