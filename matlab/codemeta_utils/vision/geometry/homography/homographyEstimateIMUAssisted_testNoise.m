function homographyEstimateIMUAssisted_testNoise
resetRands(1)
[X,G0,NVec,idxX]=homFlowDatasetStructure('NPlanes',2);
[GWorldTruth,xTruth]=homFlowDatasetMotion(X,G0,'NNewCameras',3);

GTruth=computeRelativePoseFromG(GWorldTruth(:,:,1),GWorldTruth(:,:,2:end),'references','12');
[RTruth,TTruth]=G2RT(GTruth);
NVec1Truth=rigidTransformG(GWorldTruth(:,:,1),NVec,'references','wc','planes');
[N1Truth,d1Truth]=planeNVecToNd(NVec1Truth);

%optsEstimation={};
optsEstimation={'flagCompensateRotation',true,'maxIt',2};

sigmaRotation=(0:5)*0.01;
sigmaImages=(0:3)*0.01;

NTrials=30;
NSigmaRotation=length(sigmaRotation);
NSigmaImages=length(sigmaImages);

errors.T=zeros(NSigmaRotation,NSigmaImages,NTrials);
errors.N=zeros(NSigmaRotation,NSigmaImages,NTrials);
for iTrial=1:NTrials
    fprintf('# Trial %d/%d\n',iTrial,NTrials)
    noiseRotations=randn(3,size(RTruth,3));
    noiseImages=randn(size(xTruth));
    for iSigmaRotation=1:NSigmaRotation
        fprintf('## SigmaRotation %d/%d\n',iSigmaRotation,NSigmaRotation)
        RRel=rot_exp(RTruth,rot_hat(RTruth,sigmaRotation(iSigmaRotation)*noiseRotations));
        for iSigmaImages=1:NSigmaImages
            fprintf('### SigmaImages %d/%d\n',iSigmaImages,NSigmaImages)
            x=xTruth+sigmaImages(iSigmaImages)*noiseImages;

            [TEst,NVec1Est]=homographyEstimateIMUAssisted(x,RRel,idxX,optsEstimation{:});
            [N1Est,d1Est]=planeNVecToNd(NVec1Est);
            s=d1Est(1)/d1Truth(1);

            errors.N(iSigmaRotation,iSigmaImages,iTrial)=mean(sqrt(sum((N1Truth-N1Est).^2)));
            errors.T(iSigmaRotation,iSigmaImages,iTrial)=mean(sqrt(sum((TEst/s-TTruth).^2)));
        end
    end
end

fileNameSave=mfilename;
if ~isempty(optsEstimation)
    fileNameSave=[fileNameSave '_' cell2concat(optsEstimation)];
end
fileNameSave=[fileNameSave datestr(now)];
fileNameSave=fileNameClean(fileNameSave);
save(fileNameSave)
disp(['Results saved to ' fileNameSave])

