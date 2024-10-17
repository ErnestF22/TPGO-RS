function sfm_rawEssentialRansac_test
NOutliers=30;
[G1,G2,x1,x2]=epipolar_dataset();
ETruth=epipolarBuildEFromG(G1,G2,'references');
x1=[rand(2,NOutliers) x1];
x2=[rand(2,NOutliers) x2];
[EEst,~,bestflagInlier] = sfm_rawEssentialRansac(x1,x2,100,1e-5,'methodE','5pt+validation');
disp([ETruth EEst])
disp(bestflagInlier)
