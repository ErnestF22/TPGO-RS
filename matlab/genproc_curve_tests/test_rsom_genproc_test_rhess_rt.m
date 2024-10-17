function test_rsom_genproc_test_rhess_rt
problem=test_rsom_genproc();
% curveR=test_rsom_genproc_curve(problem);

T_rand = rand(3,5);
vT_rand = rand(3,5);
[c,dc,~,~,ddc]=real_geodFun(T_rand, vT_rand);

f=@(t) problem.gs2(c(t));
gradf=@(t) problem.htr(dc(t));
df=@(t) sum(stiefel_metric(c,gradf(t),dc(t)));
funCheckDer(f,gradf)

end %file function

