%function x=quadlinsys(z,Q,W,tol)
%Solves for x in linear systems of the kind Q*W*Q'*x=z where Q has more row
%than columns and is possibly rank-deficient
%If omitted, W is the identity matrix of opportune dimension
%If present, tol specify the tolerance on the singular values of Q to
%compute its rank
function x=quadlinsys(z,Q,W,tol)
if(exist('W','var')==0)
    W=eye(size(Q,2));
end
if(exist('tol','var')==0)
    tol=0;
end

[U,S,V]=svd(Q');

r=sum(diag(S)>tol);

U=V(:,1:r);
S=S(1:r,1:r);
V=U(:,1:r);

x=(S*V'*W*V*S*U')\(U'*z);
