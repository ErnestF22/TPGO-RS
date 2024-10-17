function sphere3_logDiffMat_test
[e1,de1]=sphere_randGeodFun();
[e2,de2]=sphere_randGeodFun();
de=@(t) [de1(t);de2(t)];

l=@(t) sphere_log(e1(t),e2(t));
Dl=@(t) sphere3_logDiffMat(e1(t),e2(t));
dl=@(t) Dl(t)*de(t);

funCheckDer(l,dl,'angle')
