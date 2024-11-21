function test_check_ssom_test_egrad()
problem=test_check_ssom();
% curve=test_check_ssom_curve(problem);

% N = problem.sz(3);
% nrs = problem.sz(1);
% d = problem.sz(2);
e = size(problem.edges,1);

lambda0=rand(e,1,1);
vLambda0=rand(e,1,1);

[lambda,dLambda,~,~,~]=real_geodFun(lambda0, vLambda0);
curve.c=@(t) lambda(t);
curve.dc=@(t) dLambda(t);

f=@(t) problem.cost(curve.c(t));
egradf=@(t) problem.egrad_lambda(curve.c(t));
df=@(t) sum(stiefel_metric([],egradf(t),curve.dc(t)));
funCheckDer(f,df)