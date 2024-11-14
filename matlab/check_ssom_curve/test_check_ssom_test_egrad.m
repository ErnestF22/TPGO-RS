function test_check_ssom_test_egrad()
problem=test_check_ssom();
curve=test_check_ssom_curve(problem);

sz=problem.sz;
% sz(3) = 1; %remove this later!
R10=rand(8,1,1);
vR10=rand(8,1,1);

[R1,dR1,~,~,ddR1]=real_geodFun(R10, vR10);
curve.c=@(t) R1(t);
curve.dc=@(t) dR1(t);
% funCheckDer(multitrace(curve.c), multitrace(curve.dc))
curve.ddc=@(t) ddR1(t);

f=@(t) problem.cost(curve.c(t));
egradf=@(t) problem.egrad_lambda(curve.c(t));
df=@(t) sum(stiefel_metric(curve.c,egradf(t),curve.dc(t)));
funCheckDer(f,df)