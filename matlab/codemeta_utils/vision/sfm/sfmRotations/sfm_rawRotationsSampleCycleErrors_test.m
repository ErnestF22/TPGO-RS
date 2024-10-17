function sfm_rawRotationsSampleCycleErrors_test
%load sfm_test_data_fountain
load sfm_test_data_synthetic_clean.mat
E=sfm_matchEdges(data,'memberMatch','matchFiltered');
w=sfm_matchCounts(data,'memberMatch','matchFiltered');
w=max(w)-w;
R=G2R(data.matchPoseEstimated);
[C,e]=sfm_rawRotationsSampleCycleErrors(R,[E w]);
l=sum(abs(C));
[lUnique,eMean]=funIndexed(l,e,@mean);
plot(lUnique,eMean*180/pi)
