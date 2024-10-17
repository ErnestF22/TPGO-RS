function homographyEstimateIMUAssisted_test

resetRands(2)

%generate structure and motion
[X,G0,NVec,idxX]=homFlowDatasetStructure('NPlanes',2);
[GWorldTruth,xTruth]=homFlowDatasetMotion(X,G0,'NNewCameras',3);

%show ground truth in 3-D and the generated images
figure(1)
subplot(1,2,1)
homFlowDatasetShow(X,GWorldTruth,NVec,idxX)
subplot(1,2,2)
homFlowDatasetShow(X,GWorldTruth,NVec,idxX)
view(-90,0)
%keyboard
figure(2)
homFlowDatasetShowImages(xTruth,idxX)

%generate ground-truth parameters for the homography
sigmaRotation=0.03;
sigmaImages=0.00;
GTruth=computeRelativePoseFromG(GWorldTruth(:,:,1),GWorldTruth(:,:,2:end),'references','12');
[RTruth,TTruth]=G2RT(GTruth);
NVec1Truth=rigidTransformG(GWorldTruth(:,:,1),NVec,'references','wc','planes');
[N1Truth,d1Truth]=planeNVecToNd(NVec1Truth);

%perturbe rotations and images
RRel=rot_randn(RTruth,sigmaRotation);
x=xTruth+sigmaImages*randn(size(xTruth));

%run algorithm
[TEst,NVec1Est,RRelEst]=homographyEstimateIMUAssisted(x,RRel,idxX,...
    'flagCompensateRotation',true,'maxIt',2);

%extract estimated plane in normal and distance
[N1Est,d1Est]=planeNVecToNd(NVec1Est);
%scale for comparison purposes
s=d1Est(1)/d1Truth(1);

disp('[N1Truth N1Est]')
disp([N1Truth N1Est])
disp(meanError(N1Truth,N1Est));
disp('[TEst/s TTruth]');
disp([TEst/s TTruth]);
disp(meanError(TEst/s,TTruth));

disp('Mean rotation error')
disp(mean(rot_dist(RTruth,RRelEst,'vector')))


function e=meanError(A,B)
e=mean(sqrt(sum((A-B).^2)));