function test_problem_test_egrad()
problem=test_problem();
curve=test_problem_curve(problem);

c=curve.c;
dc=curve.dc;

f=@(t) problem.cost(c(t));
egradf=@(t) problem.egrad(c(t));
df=@(t) sum(stiefel_metric(c,egradf(t),dc(t)));
funCheckDer(f,df)