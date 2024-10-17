%function PvOrth=orthComplementProjector(v)
%Compute a projection matrix onto the orthogonal complement of v, i.e.,
%PvOrth*v=0.
function PvOrth=orthComplementProjector(v)
[U,S,~]=svd(v);
if min(size(S))==1
    s=S(1,1);
else
    s=diag(S);
end
tol = max(size(v)) * eps(max(s));
r=sum(s>tol);

PvOrth=U(:,r+1:end)*U(:,r+1:end)';
