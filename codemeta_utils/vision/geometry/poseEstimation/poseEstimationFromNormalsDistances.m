%Estimate pose using plane (normals and distances) correspondences
%function [R,T]=poseEstimationFromNormalsDistances(n1,d1,n2,d2)
%The planes are supposed to have equations n'*x=d. The function returns the
%rigid transformation mapping points in reference 2 to reference 1.
function [R,T]=poseEstimationFromNormalsDistances(n1,d1,n2,d2)
[U,~,V]=svd(n1*n2');
R=U*diag([1,1,det(U*V')])*V'; %use determinant to fix minimal case with only two pairs
if size(n1,2)==3
    T=pinv(n1')*(d1-d2)';
else
    T=n1'\(d1-d2)';
end
