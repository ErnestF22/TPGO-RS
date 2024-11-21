function test_check_ssom_test_ehess_lambda_T()
problem=test_check_ssom();
curve=test_check_ssom_curve(problem);

c=curve.c;
dc=curve.dc;

% f=@(t) problem.cost(c(t));
% gradf=@(t) problem.grad(c(t));
% df=@(t) sum(stiefel_metric(c,gradf(t),dc(t)));
% funCheckDer(f,df)

