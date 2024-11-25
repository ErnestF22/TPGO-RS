function test_check_ssom_test_ehess_T_lambda()
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
[T,dT,~,~,ddT]=real_geodFun(T0, vT0);
[R,~,~,~,~]=rot_geodFun(R0, vR0);

curve.c=@(t) T(t);
curve.dc=@(t) dT(t);
curve.ddc=@(t) ddT(t);

% f=@(t) problem.cost(curve.c(t));
egradf=@(t) problem.egrad_T(curve.c(t),R(t),lambda(t));
df=@(t) sum(stiefel_metric([],egradf(t),curve.dc(t)));
% funCheckDer(f,df)
ehessf = @(t) problem.ssom_ehess_t_lambda(curve.c(t),curve.dc(t),R(t), dLambda(t));
ddf_1 = @(t) stiefel_metric([], ehessf(t), curve.dc(t), 'euclidean');
ddf_2 = @(t) stiefel_metric([], egradf(t), curve.ddc(t), 'euclidean');
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



