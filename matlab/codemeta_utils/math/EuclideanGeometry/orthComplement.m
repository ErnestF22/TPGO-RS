%function vOrth=orthComplement(v)
%Compute a basis for the orthogonal complement to the columns of v, i.e.,
%vOrth'*v=0
function vOrth=orthComplement(v)
r=size(v,2);
[U,S,V]=svd(v);
vOrth=U(:,r+1:end);
