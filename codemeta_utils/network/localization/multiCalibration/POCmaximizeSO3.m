function POCmaximizeSO3
resetRands()
R0=rot_randn();
T0=randn(3,1);
NPlanes=5;
n1=cnormalize(randn(3,NPlanes));
d1=rand(1,NPlanes);
n2=permute(sphere_randn(permute(R0'*n1,[1 3 2]),0.1),[1 3 2]);
d2=d1-T0'*n1;

R=poseEstimationFromNormalsDistances(n1,d1,n2,d2);
disp([R R0])

[eR,~,gradRVec]=poseEstimationFromNormalsDistancesResiduals(R,T0,n1,d1,n2,d2,'methodR','cosine');
disp(sum(1-eR))
disp(sum(gradRVec,2))
