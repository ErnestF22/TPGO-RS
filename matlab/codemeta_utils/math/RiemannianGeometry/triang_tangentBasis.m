function T=triang_tangentBasis(x)
 vOrth=triang_tangentBasisOrth(x);
 O=orthCompleteBasis(vOrth);
 T=O(:,2:end);
 
 