function test_check_ssom_test_egrad()


problem=test_check_ssom();
% curve=test_check_ssom_curve(problem);

N = problem.sz(3);
nrs = problem.sz(1);
d = problem.sz(2);
e = size(problem.edges,1);

lambda0=rand(e,1,1);
vLambda0=rand(e,1,1);
[lambda,dLambda,~,~,~]=real_geodFun(lambda0, vLambda0);

R0=randrot(d,N);
vR0=rot_randTangentNormVector(R0);
[R,dR,~,~,~,~]=rot_geodFun(R0, vR0);

T0 = rand(nrs,N);
vT0 = rand(nrs,N);
[T,dT,~,~,~]=real_geodFun(T0, vT0);

%% lambda

disp("grad lambda")
curve.c=@(t) lambda(t);
curve.dc=@(t) dLambda(t);
f=@(t) problem.cost(R0, T0, curve.c(t));
egradf=@(t) problem.egrad_lambda(R0,T0,curve.c(t));
df=@(t) sum(stiefel_metric([],egradf(t),curve.dc(t)));
funCheckDer(f,df)

%% T
disp("grad T")

curve.c=@(t) T(t);
curve.dc=@(t) dT(t);

f=@(t) problem.cost(R0, curve.c(t), lambda0);
egradf=@(t) problem.egrad_T(R0,curve.c(t),lambda0);
df=@(t) sum(stiefel_metric([],egradf(t),curve.dc(t)));
funCheckDer(f,df)

%% R

disp("grad R")

curve.c=@(t) R(t);
curve.dc=@(t) dR(t);

f=@(t) problem.cost(curve.c(t), T0, lambda0);
egradf=@(t) problem.egrad_R(curve.c(t), T0, lambda0);
df=@(t) sum(stiefel_metric([],egradf(t),curve.dc(t)));
% funCheckDer(f,df)
funCheckDer(f,df)

end %file function

