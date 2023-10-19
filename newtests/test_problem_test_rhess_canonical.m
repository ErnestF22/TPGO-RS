function test_problem_test_rhess_canonical()
problem=test_problem();
curve=test_problem_curve(problem);

c=curve.c;
dc=curve.dc;
ddc=curve.ddc;

f=@(t) problem.cost(c(t));
rgradf=@(t) problem.rgrad(c(t));
df=@(t) sum(stiefel_metric(c(t),rgradf(t),dc(t),'canonical'));

%%%
rhessf = @(t) problem.rhess(c(t));
ddf_1 = @(t) sum(stiefel_metric(c(t), rhessf(t), dc(t), 'canonical'));
ddot_c = @(t) stiefel_tangentProj(c(t), ddc(t));
ddf_2 = @(t) sum(stiefel_metric(c(t), rgradf(t), ddot_c(t), 'canonical'));
ddf = @(t) ddf_1(t) + ddf_2(t);

figure('Name','Riem. Hessian: cost along curve')
funCheckDer(df,ddf,linspace(-1,1))