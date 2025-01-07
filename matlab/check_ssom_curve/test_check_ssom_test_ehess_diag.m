function test_check_ssom_test_ehess_diag()
%%

problem=test_check_ssom();
% curve=test_check_ssom_curve(problem);

N = problem.sz(3);
% nrs = problem.sz(1);
d = problem.sz(2);
e = size(problem.edges,1);

lambda0=rand(e,1,1);
vLambda0=rand(e,1,1);

R0 = eye3d(d,d,N);
vR0 = eye3d(d,d,N);

[lambda,dLambda,~,~,ddLambda]=real_geodFun(lambda0, vLambda0);
[R,dR,~,~,~,ddR]=rot_geodFun(R0, vR0);
T0 = rand(d,N);
vT0 = rand(d,N);
[T,dT,~,~,ddT]=real_geodFun(T0, vT0);


%% lambda lambda
curve.c=@(t) lambda(t);
curve.dc=@(t) dLambda(t);
curve.ddc=@(t) ddLambda(t);

% f=@(t) problem.cost(curve.c(t));
egradf=@(t) problem.grad_lambda(curve.c(t),R(t),T0);
df=@(t) sum(stiefel_metric([],egradf(t),curve.dc(t)));
% funCheckDer(f,df)
ehessf = @(t) problem.ssom_ehess_lambda_lambda(curve.c(t),curve.dc(t),R(t));
ddf_1 = @(t) stiefel_metric([], ehessf(t), curve.dc(t), 'euclidean');
ddf_2 = @(t) stiefel_metric([], egradf(t), curve.ddc(t), 'euclidean');
ddf = @(t) sum(ddf_1(t) + ddf_2(t));    

funCheckDer(df,ddf)

%% T T
curve.c=@(t) T(t);
curve.dc=@(t) dT(t);
curve.ddc=@(t) ddT(t);

% f=@(t) problem.cost(curve.c(t));
egradf=@(t) problem.egrad_T(curve.c(t));
df=@(t) sum(stiefel_metric([],egradf(t),curve.dc(t)));
% funCheckDer(f,df)
ehessf = @(t) problem.ssom_ehess_T_T(curve.c(t),curve.dc(t));
ddf_1 = @(t) stiefel_metric([], ehessf(t), curve.dc(t), 'euclidean');
ddf_2 = @(t) stiefel_metric([], egradf(t), curve.ddc(t), 'euclidean');
ddf = @(t) sum(ddf_1(t) + ddf_2(t));    

funCheckDer(df,ddf)

%% T T
curve.c=@(t) R(t);
curve.dc=@(t) dR(t);
curve.ddc=@(t) ddR(t);

% f=@(t) problem.cost(curve.c(t));
egradf=@(t) problem.rgrad_R(curve.c(t));
df=@(t) sum(stiefel_metric([],egradf(t),curve.dc(t)));
% funCheckDer(f,df)
ehessf = @(t) problem.ssom_rhess_R_R(curve.c(t),curve.dc(t));
ddf_1 = @(t) stiefel_metric([], ehessf(t), curve.dc(t), 'euclidean');
ddf_2 = @(t) stiefel_metric([], egradf(t), curve.ddc(t), 'euclidean');
ddf = @(t) sum(ddf_1(t) + ddf_2(t));    

funCheckDer(df,ddf)



end %file function
