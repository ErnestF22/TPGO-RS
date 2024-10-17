%Refine pose estimation from poseEstimationN
%function [R,T]=poseEstimationFromNormalsPointsRefine(n,d,x2D,R)
%Minimizes the squared norm of the residuals
function [R,T]=poseEstimationFromNormalsPointsRefine(n,d,x2D,R)
M=n';
[UM,~,~]=svd(M,'econ');
NPoints=size(M,1);
nx=zeros(6,NPoints);
for inx=1:NPoints
    nx(:,inx)=reshape(n(:,inx)*x2D(:,inx)',6,1);
end
p=d'-UM*(UM'*d');
A=[zeros(NPoints,3) nx'-UM*(UM'*nx')];
R=rot_minimizeLinearFunctional(A,p,R);
T=poseEstimationFromNormalsPoints_TFromR(R,nx,M,d);
