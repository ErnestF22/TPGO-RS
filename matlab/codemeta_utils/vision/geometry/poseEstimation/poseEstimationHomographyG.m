function [G21,N1,lambda]=poseEstimationHomographyG(x1,x2)
HEst=homographyNormalize(homographyEstimateLinear(x1,x2));
[G21,N1,lambda]=homographyToG(HEst,x1,x2);
