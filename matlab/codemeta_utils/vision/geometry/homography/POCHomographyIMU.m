function POCHomographyIMU
%resetRands(1)
[X,G0,NVec,idxX]=homFlowDatasetStructure('NPlanes',2);
[GWorldTruth,xTruth]=homFlowDatasetMotion(X,G0,'NNewCameras',3);

figure(1)
homFlowDatasetShow(X,GWorldTruth,NVec,idxX)
figure(2)
homFlowDatasetShowImages(xTruth,idxX)

sigmaRotation=0.00;
sigmaImages=0.00;
GTruth=computeRelativePoseFromG(GWorldTruth(:,:,1),GWorldTruth(:,:,2:end),'references','12');
[RTruth,TTruth]=G2RT(GTruth);
NVec1Truth=rigidTransformG(GWorldTruth(:,:,1),NVec,'references','wc','planes');
[N1Truth,d1Truth]=planeNVecToNd(NVec1Truth);

RRel=rot_randn(RTruth,sigmaRotation);
x=xTruth+sigmaImages*randn(size(xTruth));

[TEst,NVec1Est,w]=homographyEstimateIMUAssisted(x,idxX,RRel);

[N1Est,d1Est]=planeNVecToNd(NVec1Est);

s=d1Est(1)/d1Truth(1);

disp([N1Truth N1Est])
disp([TEst/s TTruth]);

disp(w)

function [T,NVec,w]=homographyEstimateIMUAssisted(x,idxX,RRel)
NFrames=size(x,3)-1;
NPlanes=max(idxX);
NPairs=NFrames*NPlanes;
HContinuous=zeros(3,3,NPairs);
E=zeros(NPairs,2);
idxPlanes=reshape(1:NPairs,NFrames,NPlanes);
for iPlane=1:NPlanes
    xGroup=x(:,idxX==iPlane,:);
    H=homographyEstimateLinear(xGroup);
    H=homographyNormalize(H,xGroup);
    HContinuous(:,:,idxPlanes(:,iPlane))=multiprod(multitransp(RRel),H)-repmat(eye(3),1,1,size(H,3));
    E(idxPlanes(:,iPlane),1)=1:NFrames;
    E(idxPlanes(:,iPlane),2)=iPlane;
end
[w,v,n]=homographyContinuousEstimateNuclearNorm(HContinuous,E);
T=-squeeze(multiprod(RRel,permute(v,[1 3 2])));
NVec=planeNScaledToNVec(n);

