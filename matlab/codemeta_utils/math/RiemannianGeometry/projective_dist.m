%Compute distance in projective space
%function D=projective_dist(X,Y)
%The distance is based on the Riemannian distance in the sphere
function D=projective_dist(X,Y)
Dp=sphere_dist(X,Y);
Dn=sphere_dist(X,-Y);
D=min(Dn,Dp);
