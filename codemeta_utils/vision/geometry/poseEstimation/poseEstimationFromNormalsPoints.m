%Estimate pose from planes (normals and distances) to 3-D points in a common plane
%function [R,T]=poseEstimationFromNormalsPoints(n,d,x2D)
%The planes are supposed to have equations n1'*x=d1 in the reference frame 1.
%The points are supposed to be of the form x=[zeros(); x2D], i.e., only the
%y and z coordinates are passed as inputs, and the x coordinates are
%assumed to be zero
%The function returns the rigid transformation mapping points in reference
%2 to reference 1. At least 6 independent correspondences are necessary
function [R,T]=poseEstimationFromNormalsPoints(n,d,x2D)
M=n';
[UM,~,~]=svd(M,'econ');
NPoints=size(M,1);
nx=zeros(6,NPoints);
for inx=1:NPoints
    nx(:,inx)=reshape(n(:,inx)*x2D(:,inx)',6,1);
end
p=d'-UM*(UM'*d');
A=nx'-UM*(UM'*nx');
RPartEst=reshape(A\p,3,2);
[UR,~,VR]=svd([zeros(3,1) RPartEst]);
R=UR*diag([1 1 det(UR*VR')])*VR';
T=poseEstimationFromNormalsPoints_TFromR(R,nx,M,d);
