%Compute errors on plane (normals and distances) correspondences
function [eRot,eTransl]=poseEstimationFromNormalsDistancesResidualsG(G,n1,d1,n2,d2)
[R,T]=G2RT(G);
[eRot,eTransl]=poseEstimationFromNormalsDistancesResiduals(R,T,n1,d1,n2,d2);
