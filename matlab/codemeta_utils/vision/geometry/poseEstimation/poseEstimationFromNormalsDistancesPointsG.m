%Estimate pose from planes (normals and distance) to 3-D point-pairs correspondences
%function [R,T]=poseEstimationFromNormalsDistancesPointsG(n1,d1,x2)
%Same as poseEstimationFromNormalsDistancesPoints, but returns the single
%pose G instead of the pair R,T
function G=poseEstimationFromNormalsDistancesPointsG(n1,d1,x2)
[R,T]=poseEstimationFromNormalsDistancesPoints(n1,d1,x2);
G=RT2G(R,T);