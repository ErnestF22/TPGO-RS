%Compute derivative of bearing vectors
%function [dy,Q]=bearingComputeDerivative(dx,y,ny)
%Inputs
%   dx   derivative of robot's location
%   y   current bearings
%   ny  ranges corresponding to the bearings
%Output
%   dy  derivative of Y along v
%   Q   array of matrices such that dY=Q*v
function [dy,Q]=bearingComputeDerivative(dx,y,ny)
flagComputeQ=false;
if nargout>1
    flagComputeQ=true;
end
[d,Ny]=size(y);
dy=zeros(d,Ny);
if flagComputeQ
    Q=zeros(d,d,Ny);
end
I=eye(d);
for iY=1:Ny
    Qi=-inv(ny(iY))*(I-y(:,iY)*y(:,iY)');
    dy(:,iY)=Qi*dx;
    if flagComputeQ
        Q(:,:,iY)=Qi;
    end
end
