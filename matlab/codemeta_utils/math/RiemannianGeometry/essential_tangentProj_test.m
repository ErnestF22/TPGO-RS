function essential_tangentProj_test
[Qt,vt]=essential_randGeodFun();

pt=@(t) essential_tangentProjVerticalScalar(Qt(t),vt(t));
 
plotfun(pt,linspace(0,4*pi))
