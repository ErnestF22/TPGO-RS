%Compute second derivative of bearing vectors
%function ddy=bearingComputeSecondDerivative(dx,ddx,y,ny)
%Inputs
%   dx  derivative of robot's location (i.e., \dot{x})
%   ddx second derivative of location
%   Y   current bearings
%   nY  ranges corresponding to the bearings
%Output
%   ddY  second derivative of Y along v
function ddy=bearingComputeSecondDerivative(dx,ddx,y,ny)
% flagComputeQ=false;
% if nargout>1
%     flagComputeQ=true;
% end
[d,Ny]=size(y);
dy=bearingComputeDerivative(dx,y,ny);
ddy=zeros(d,Ny);
% if flagComputeQ
%     Q=zeros(d,d,NY);
% end
I=eye(d);
for iY=1:Ny
    yi=y(:,iY);
    dyi=dy(:,iY);
    ri=ny(iY);
    Py=I-yi*yi';
    ddy(:,iY)=-ri^-2*dx'*yi*Py*dx+ri^-1*(dyi*yi'+yi*dyi')*dx-ri^-1*Py*ddx;
%     if flagComputeQ
%         Q(:,:,iY)=Qi;
%     end
end
