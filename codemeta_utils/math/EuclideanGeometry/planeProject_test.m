function planeProject_test
X=randn(3,1000);
NVec=planeRandn(eye(4,1));
XProj=planeProject(NVec,X);
disp(norm(planeResiduals(NVec,XProj),Inf))

