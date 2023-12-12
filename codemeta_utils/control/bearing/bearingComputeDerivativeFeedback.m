function [dx,Q] = bearingComputeDerivativeFeedback(dy,y,ny)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
[d,Ny]=size(y);
I=eye(d);
for iY=1:Ny
    Qi=-inv(ny(iY))*(I-y(:,iY)*y(:,iY)');
    Q(1:d,(iY-1)*d+1:d*iY)=Qi;
end
dy=reshape(dy,12,1);
dx=Q*dy;

