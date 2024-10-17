function Sigma=poseEstimationCovarianceFromG(G,X,x)
Sigma=poseEstimationCovariance(G2R(G),G2T(G),X,x);