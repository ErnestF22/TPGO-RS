function test_rsom_test_rgrad()
problem=test_rsom();
curve=test_rsom_curve(problem);

c=curve.c;
dc=curve.dc;

f=@(t) problem.cost(c(t));
gradf=@(t) problem.rgrad(c(t));
df=@(t) sum(stiefel_metric(c,gradf(t),dc(t)));
funCheckDer(f,df)

