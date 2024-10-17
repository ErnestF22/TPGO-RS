%Estimate pose using plane (normals and distances) correspondences
%function G=poseEstimationFromNormalsDistances(n1,d1,n2,d2)
%Same as poseEstimationFromNormalsDistances, but returns pose G instead of
%R and T separately.
function G=poseEstimationFromNormalsDistancesG(n1,d1,n2,d2)
[R,T]=poseEstimationFromNormalsDistances(n1,d1,n2,d2);
G=RT2G(R,T);
