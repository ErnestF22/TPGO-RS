%Estimate pose from planes (normals and distances) to 3-D points in a common plane
%function G=poseEstimationFromNormalsPointsG(n,d,x2D)
%Same as poseEstimationFromNormalsPoints but returns a single pose G
function G=poseEstimationFromNormalsPointsG(n,d,x2D)
[R,T]=poseEstimationFromNormalsPoints(n,d,x2D);
G=RT2G(R,T);
