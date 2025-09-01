function test_check_ssom_test_ehess_R_T()

problem=test_check_ssom();
% curve=test_check_ssom_curve(problem);

N = problem.sz(3);
% nrs = problem.sz(1);
d = problem.sz(2);
e = size(problem.edges,1);

lambda0=rand(e,1,1);
vLambda0=rand(e,1,1);

R0 = make_rand_stiefel_3d_array(d,d,N); %first d with nrs later
% vR0 = eye3d(d,d,N);

T0 = rand(d,N);
vT0 = rand(d,N);

% [lambda,dLambda,~,~,ddLambda]=real_geodFun(lambda0, vLambda0);
[T,dT,~,~,ddT]=real_randGeodFun(T0);
% [R,dR,~,~,ddR]=rot_geodFun(R0, vR0);
% 
curve.c=@(t) T(t);
curve.dc=@(t) dT(t);
curve.ddc=@(t) ddT(t);

% f=@(t) problem.cost(curve.c(t));
gradf=@(t) vec(problem.rgrad_R(R0, curve.c(t), lambda0));
% df=@(t) sum(stiefel_metric([],gradf(t),curve.dc(t)));
% funCheckDer(f,df)
ehessf = @(t) vec(problem.ssom_ehess_R_T(R0, curve.c(t), curve.dc(t), lambda0));
% ddf_1 = @(t) stiefel_metric([], ehessf(t), curve.dc(t), 'euclidean');
% ddf_2 = @(t) stiefel_metric([], gradf(t), curve.ddc(t), 'euclidean');
ddf = @(t) sum(ddf_1(t) + ddf_2(t));    

% funCheckDer(df,ddf,'angle')
funCheckDer(gradf, ehessf)

end %file function