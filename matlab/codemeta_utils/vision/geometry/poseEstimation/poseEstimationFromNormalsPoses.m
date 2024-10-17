%Estimate pose from pose-plane correspondences
%function [R,T]=poseEstimationFromNormalsPoses(n,d,G)
%The pair n,d describes the plane in the reference frame (reference
%interpretation) which we need to estimate, while G represents the pose of
%a reference frame 
function [R,T]=poseEstimationFromNormalsPoses(n,d,Rn,Tn)
NPlanes=size(n,2);
C=zeros(3,3);
A=zeros(NPlanes,3);
b=zeros(NPlanes,1);
for iPlane=1:NPlanes
    r=Rn(:,3,iPlane);
    %for the rotations
    C=C+r*n(:,iPlane)';
    %for the translations
    A(iPlane,:)=r';
    b(iPlane)=r'*Tn(:,iPlane)-d(iPlane);
end
[U,S,V]=svd(C);
R=U*diag([1,1,det(U*V')])*V';
T=A\b;