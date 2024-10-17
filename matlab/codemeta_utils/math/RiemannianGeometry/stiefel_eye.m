%returns the identity element for the Stiefel manifold
function y=stiefel_eye(y)
[n,p]=size(y);
y=eye(n);
y=y(:,1:p);
