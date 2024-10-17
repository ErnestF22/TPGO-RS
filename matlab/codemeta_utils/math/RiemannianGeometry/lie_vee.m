%function vvec=lie_vee(lf,y,v)
%Given vectors v in the tangent space at y, gives their cooridnate
%representation in a base for the tangent space
function vvec=lie_vee(lf,y,v)
T=lf.tangentBasis(y);
vvec=lf.metric(y,T,v);
