function sfm_rawRotationsCycleErrors2EdgeLabels_test
%resetRands()

%dataset
load sfm_test_data_synthetic_clean.mat
E=sfm_matchEdges(data,'memberMatch','matchFiltered');
w=sfm_matchCounts(data,'memberMatch','matchFiltered');
w=max(w)-w;
NEdges=size(E,1);
R=G2R(data.matchPoseEstimated);
xeTruth=zeros(NEdges,1);

%generate outliers
NOutliers=3;
v=[45 90]*pi/180;
idxOutliers=randperm(NEdges,NOutliers);
xeTruth(idxOutliers)=1;
R(:,:,idxOutliers)=rot_randNotch(R(:,:,idxOutliers),v,[],'U',diag([1 1 0]));

%sample cycles
[C,e]=sfm_rawRotationsSampleCycleErrors(R,[E w],'triangles','NSamples',100);
NCycles=size(C,2);
xlTruth=zeros(NCycles,1);
for iCycles=1:NCycles
    c=C(:,iCycles)>0;
    xlTruth(iCycles)=max(xeTruth(c));
end

%do inference
xe=sfm_rawRotationsCycleErrors2EdgeLabels(C,e);

disp([xe xeTruth])
err=sum(xe~=xeTruth);
fprintf('%d errors\n',err)

