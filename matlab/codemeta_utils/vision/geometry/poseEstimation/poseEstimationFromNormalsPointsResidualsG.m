function e=poseEstimationFromNormalsPointsResidualsG(G,n,d,x2D)
[R,T]=G2RT(G);
e=poseEstimationFromNormalsPointsResiduals(R,T,n,d,x2D);
