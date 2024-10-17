function v=triang_tangentProject(x,v)
vOrth=triang_tangentBasisOrth(x);
v=v(:)-(v(:)'*vOrth)*vOrth;
v=reshape(v,size(x));
