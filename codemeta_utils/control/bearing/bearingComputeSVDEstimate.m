%Estimate dx and distances to landmarks for dy and Q matrices.
% function [dx,ny]=bearingComputeSVDEstimate(y,dy)
%Inputs
%   y current bearings
%   dy current derivative of bearings (bearingComputeDerivative)
%Outputs
%   dx current estimated derivative of robot's location
%   ny current estimated distances from landmarks
function [landmarks] = bearingComputeSVDEstimate(y,dy,dx)
[d,Ny]=size(y);
A=zeros(d*Ny,d+Ny);
b=zeros(d+Ny,1);
I=eye(d);
for iY=1:Ny
    Qi=(I-y(:,iY)*y(:,iY)');
    A(1+d*(iY-1):d*iY,1:d)=Qi;
    A(1+d*(iY-1):d*iY,iY+d)=dy(:,iY);
end
[~,~,V]=svd(A,'econ');
NS=V(:,Ny+1:end);
alpha=NS(1:d,:)\dx;
v1=NS*alpha;
landmarks=y.*v1(d+1:end)';


