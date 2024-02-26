function test_rsom_genproc_test_rhess_tr
problem=test_rsom_genproc();
% curveR=test_rsom_genproc_curve(problem);

R_rand = matUnstack([ ...
    [randrot(3)]; [randrot(3)]; [randrot(3)]; [randrot(3)]; [randrot(3)]]);
vR_rand = matUnstack([ ...
    [randrot(3)]; [randrot(3)]; [randrot(3)]; [randrot(3)]; [randrot(3)]]);
[c,dc,R0,dR0,vVec,ddc,dvVec]=rot_geodFun(R_rand, vR_rand);


f=@(t) problem.gs2(c(t));
gradf=@(t) problem.htr(dc(t));
df=@(t) sum(stiefel_metric(c,gradf(t),dc(t)));
funCheckDer(f,gradf)

end %file function
