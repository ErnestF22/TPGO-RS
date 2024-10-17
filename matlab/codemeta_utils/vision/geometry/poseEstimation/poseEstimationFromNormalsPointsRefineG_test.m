function poseEstimationFromNormalsPointsRefineG_test
load('poseEstimationFromNormalsPointsRefineG_test_inputData.mat')
G=GHokuyoCamera;
n=normalPlanesCameraFlat;
d=distancePlanesCameraFlat;
x=XLinesHokuyoFlat2D;
x=x+0.1*randn(size(x));
GEst=poseEstimationFromNormalsPointsRefineG(n,d,x,G);
[c,gradVec]=poseEstimationFromNormalsPointsCostG(GEst,n,d,x);
disp(c)
disp(sum(gradVec,2))

