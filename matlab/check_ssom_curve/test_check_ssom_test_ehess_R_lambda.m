function test_check_ssom_test_ehess_R_lambda()

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

[lambda,dLambda,~,~,~]=real_geodFun(lambda0, vLambda0);
[T,~,~,~,~]=real_geodFun(T0, vT0);
[R,dR,~,~,ddR]=rot_geodFun(R0, vR0);

curve.c=@(t) R(t);
curve.dc=@(t) dR(t);
curve.ddc=@(t) ddR(t);

% f=@(t) problem.cost(curve.c(t));
gradf=@(t) problem.rgrad_R(curve.c(t));
df=@(t) sum(stiefel_metric([],gradf(t),curve.dc(t)));
% funCheckDer(f,df)
ehessf = @(t) problem.ssom_ehess_r_lambda(R(t),vR0,T0,lambda0);
ddf_1 = @(t) stiefel_metric([], ehessf(t), curve.dc(t), 'euclidean');
ddf_2 = @(t) stiefel_metric([], gradf(t), curve.ddc(t), 'euclidean');
ddf = @(t) sum(ddf_1(t) + ddf_2(t));    

funCheckDer(df,ddf,'angle')