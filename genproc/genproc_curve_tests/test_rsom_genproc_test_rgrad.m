function test_rsom_genproc_test_rgrad()
problem=test_rsom_genproc();
curve=test_rsom_genproc_curve(problem);

c=curve.c;
dc=curve.dc;

f=@(t) problem.cost(c(t));
gradf=@(t) problem.grad(c(t));
df=@(t) sum(stiefel_metric(c,gradf(t),dc(t)));
funCheckDer(f,df)

