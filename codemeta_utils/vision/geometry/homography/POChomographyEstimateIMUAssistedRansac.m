function POChomographyEstimateIMUAssistedRansac
resetRands(1)
sigmaRotation=0.00;
sigmaImages=0.0001;

[X,G0,NVec,idxX]=homFlowDatasetStructure('NPlanes',1);
[GWorldTruth,xTruth]=homFlowDatasetMotion(X,G0,'NNewCameras',4);

GTruth=computeRelativePoseFromG(GWorldTruth(:,:,1),GWorldTruth(:,:,2:end),'references','12');
[RTruth,TTruth]=G2RT(GTruth);
NVec1Truth=rigidTransformG(GWorldTruth(:,:,1),NVec,'references','wc','planes');
[N1Truth,d1Truth]=planeNVecToNd(NVec1Truth);


RRel=rot_randn(RTruth,sigmaRotation);
x=xTruth+sigmaImages*randn(size(xTruth));
iPlane=1;
xGroup=x(:,idxX==iPlane,:);

models=ransacModelsSample(x,@(x) funParametersEstimation(x,RRel),4,...
    'NSamples',50,'waitbar',...
    'flagUniformOutput',true);
plotPoints(models(1:3,:)); hold on; plot3(0,0,0,'o'); hold off; axis equal
hist(vec(acos(max(-1,min(1,models(1:12,:)'*models(1:12,:))))*180/pi));
xlabel('Pairwise angular distance [deg]')
save([mfilename '_data'])
%keyboard

%H=funModelEstimation(xGroup(:,1:8,:),RRel);
%plot(funResiduals(H,x))
% [H,output]=ransac(x,@(x) funModelEstimation(x,RRel),8,@funResiduals,0.01,...
%     'waitbar','Ntrials',50);
% plot(output.flagInliers)


function v=funParametersEstimation(x,RRel)
[T,NVec,R]=homographyEstimateIMUAssisted(x,RRel);
nT=norm(T(:));
T=T/nT;
NVec(4)=NVec(4)/nT;
v=[T(:);R(:);NVec];


function H=funModelEstimation(x,RRel)
[T,NVec,R]=homographyEstimateIMUAssisted(x,RRel);
H=homographyFromRT(R,T,NVec);

function e=funResiduals(H,x)
x2=homographyMap(H,x(:,:,1));
e=sqrt(sum(sum((x(:,:,2:end)-x2).^2),3))/size(x,3);

