function test_rsom_test_egrad()
problem=test_rsom();
curve=test_rsom_curve(problem);

c=curve.c;
dc=curve.dc;

f=@(t) problem.cost(c(t));
egradf=@(t) problem.egrad(c(t));
df=@(t) sum(stiefel_metric(c,egradf(t),dc(t)));
funCheckDer(f,df)