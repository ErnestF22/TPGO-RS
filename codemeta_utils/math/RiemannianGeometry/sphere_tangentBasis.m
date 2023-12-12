function T=sphere_tangentBasis(y)
d=size(y,1);
T=permute(orth(eye(d)-y*y'),[1 3 2]);
