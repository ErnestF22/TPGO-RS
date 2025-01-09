function test_check_ssom_test_ehess_lambda_R()

problem=test_check_ssom();
% curve=test_check_ssom_curve(problem);

N = problem.sz(3);
% nrs = problem.sz(1);
d = problem.sz(2);
e = size(problem.edges,1);

lambda0=rand(e,1,1);
vLambda0=rand(e,1,1);

[lambda,dLambda,~,~,ddLambda]=real_geodFun(lambda0, vLambda0);

R0=randrot(d,N);
vR0=rot_randTangentNormVector(R0);
% [R,dR,~,~,~,ddR]=rot_geodFun(R0, vR0);
[R,dR,~,~,ddR]=real_geodFun(R0, vR0);

T0 = rand(d,N);
% vT0 = rand(d,N);
% [T,dT,~,~,ddT]=real_geodFun(T0, vT0);


curve.c=@(t) R(t);
curve.dc=@(t) dR(t);
curve.ddc=@(t) ddR(t);

% f=@(t) problem.cost(curve.c(t));

grad_handle = @(t) problem.grad_lambda(lambda0, R(t), T0);
df=@(t) sum(stiefel_metric([],grad_handle(t),lambda0));
% funCheckDer(f,df)
ehessf = @(t) problem.ssom_ehess_lambda_r(dR(t), ddR(t), T0, lambda0);
ddf_1 = @(t) stiefel_metric([], ehessf(t), dLambda(t), 'euclidean');
ddf_2 = @(t) stiefel_metric([], grad_handle(t), ddLambda(t), 'euclidean');
ddf = @(t) sum(ddf_1(t) + ddf_2(t));    

funCheckDer(df,ddf)

end
