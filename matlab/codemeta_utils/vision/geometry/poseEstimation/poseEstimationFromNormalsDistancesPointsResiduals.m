%Estimate residuals on planes (normals and distance) to 3-D point-pairs correspondences
%function [eRot,eTransl]=poseEstimationFromNormalsDistancesPoints(R,T,n1,d1,x2)
%The planes are supposed to have equations n1'*x=d1 in the reference frame 1.
%The rigid transformation must map points in reference 2 to reference 1.
%Inputs
%   n1  [3 x NPlanes] normals of the planes in reference 1
%   d1  [1 x NPlanes] plane distances from the origin in reference 1
%   x2  [3 x 2 x NPlanes] pairs of 3-D points (in reference 2) belonging to
%       the planes
%Outputs
%   eRot    residuals n1'*R*l2, where l2 is the line directions extracted
%           from x2 
%   eTransl residuals n1'*(R*x2+T)-d
function [eRot,eTransl]=poseEstimationFromNormalsDistancesPointsResiduals(R,T,n1,d1,x2)
NPlanes=size(n1,2);
NPoints=2;

%convert points to line pairs
l2=zeros(3,NPlanes);
for iPlane=1:NPlanes
    l2(:,iPlane)=x2(:,:,iPlane)*[1;-1];
end
u=ones(NPoints,1);

%compute residuals
eRot=zeros(NPlanes,1);
eTransl=zeros(NPlanes,NPoints);
for iPlane=1:NPlanes
    eRot(iPlane)=n1(:,iPlane)'*R*l2(:,iPlane);
    eTransl(iPlane,:)=[1;1]*n1(:,iPlane)'*T+x2(:,:,iPlane)'*R'*n1(:,iPlane)-d1(iPlane)*u;
end
eTransl=eTransl(:);
