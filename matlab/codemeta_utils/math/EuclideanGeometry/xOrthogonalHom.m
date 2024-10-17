%Orthogonal basis to a point expressed in homogeneous coordinates
%function xOrth=xOrthogonalHom(x)
%Input
%   x       [d x N] array of points
%Output
%   xOrth   [d+1 x d x N] array of basis for the orthogonal space of
%           homogeneous(x)
function xOrth=xOrthogonalHom(x)
NX=size(x,2);
dimData=size(x,1)+1;
d=cnormalize([x;ones(1,NX)]);
d(end,:)=d(end,:)-1;
d=cnormalize(d);
I=eye(dimData,dimData-1);
xOrth=zeros(dimData,dimData-1,NX);
for iX=1:NX
    xOrth(:,:,iX)=I-2*d(:,iX)*d(1:end-1,iX)';
end

