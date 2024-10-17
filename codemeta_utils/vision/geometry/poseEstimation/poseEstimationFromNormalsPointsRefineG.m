%Refine pose estimation from poseEstimationN
%function G=poseEstimationFromNormalsPointsRefine(n,d,x2D,G)
%Same as poseEstimationFromNormalsPointsRefine, but uses a single pose
function G=poseEstimationFromNormalsPointsRefineG(n,d,x2D,G)
R=G2R(G);
[R,T]=poseEstimationFromNormalsPointsRefine(n,d,x2D,R);
G=RT2G(R,T);
