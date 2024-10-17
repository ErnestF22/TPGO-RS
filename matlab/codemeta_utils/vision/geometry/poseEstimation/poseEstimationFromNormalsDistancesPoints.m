%Estimate pose from planes (normals and distance) to 3-D point-pairs correspondences
%function [R,T]=poseEstimationFromNormalsDistancesPoints(n1,d1,x2)
%The planes are supposed to have equations n1'*x=d1 in the reference frame 1.
%The function returns the rigid transformation mapping points in reference
%2 to reference 1. At least 5 correspondences are necessary
%Inputs
%   n1  [3 x NPlanes] normals of the planes in reference 1
%   d1  [1 x NPlanes] plane distances from the origin in reference 1
%   x2  [3 x 2 x NPlanes] pairs of 3-D points (in reference 2) belonging to
%       the planes 
function [R,T]=poseEstimationFromNormalsDistancesPoints(n1,d1,x2)
NPlanes=size(n1,2);
NPoints=2;

%convert points to line pairs
l2=zeros(3,NPlanes);
for iPlane=1:NPlanes
    l2(:,iPlane)=x2(:,:,iPlane)*[1;-1];
end

%find rotation
R=rotation5pt(n1,l2);

%find translation
u=ones(NPoints,1);
idxPlanes=reshape(1:NPoints*NPlanes,NPoints,NPlanes); %indeces for each plane in the rows of M and m
M=zeros(NPoints*NPlanes,3); %M and m will contain the linear system on the translation
m=zeros(NPoints*NPlanes,1);
for iPlane=1:NPlanes
    M(idxPlanes(:,iPlane),:)=u*n1(:,iPlane)';
    m(idxPlanes(:,iPlane),:)=-x2(:,:,iPlane)'*R'*n1(:,iPlane)+d1(iPlane)*u;
end
T=M\m;
