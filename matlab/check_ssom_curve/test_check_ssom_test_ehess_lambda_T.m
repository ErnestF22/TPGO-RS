function test_check_ssom_test_ehess_lambda_T()

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
[T,~,~,~,~]=real_geodFun(T0, vT0);
[R,dR,~,~,~]=rot_geodFun(R0, vR0);

curve.c=@(t) lambda(t);
curve.dc=@(t) dLambda(t);
curve.ddc=@(t) ddLambda(t);

% f=@(t) problem.cost(curve.c(t));
egradf=@(t) problem.egrad_lambda(curve.c(t));
df=@(t) sum(stiefel_metric([],egradf(t),curve.dc(t)));
% funCheckDer(f,df)
ehessf = @(t) problem.ssom_ehess_lambda_t(curve.c(t),curve.dc(t),T(t),R(t),dR(t));
ddf_1 = @(t) stiefel_metric([], ehessf(t), curve.dc(t), 'euclidean');
ddf_2 = @(t) stiefel_metric([], egradf(t), curve.ddc(t), 'euclidean');
ddf = @(t) sum(ddf_1(t) + ddf_2(t));    

funCheckDer(df,ddf)