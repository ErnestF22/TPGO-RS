%function G=poseEstimationFromNormalsPosesG(n,d,Gn)
function G=poseEstimationFromNormalsPosesG(n,d,Gn)
[Rn,Tn]=G2RT(Gn);
[R,T]=poseEstimationFromNormalsPoses(n,d,Rn,Tn);
G=RT2G(R,T);
