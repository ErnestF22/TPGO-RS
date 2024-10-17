function c=triang_tangentProjVerticalScalar(x,v)
vOrth=triang_tangentBasisOrth(x);
c=vOrth(:)'*v(:);
