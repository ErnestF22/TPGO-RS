function vperp=lie_randTangentPerpNormVector(lf,y,v)
vperp=lf.tangentProj(y,randn(size(y)));
v=v/sqrt(lf.metric(y,v,v));
vperp=vperp-lf.metric(y,v,vperp)*v;
vperp=vperp/sqrt(lf.metric(y,vperp,vperp));
