function cameraDepthSensorPlaneCalibration_testTrials
resetRands(1)
NTrials=50;

L=3;
offset=[L/2;-L/2;0];
GTag1=RT2G(eye(3),zeros(3,1)+offset);
GTag2=RT2G(eye(3),[-L;0;0]+offset);
GTag3=RT2G(eye(3),[0;L;0]+offset);
GTag=cat(3,GTag1,GTag2,GTag3);
NIt=2000;
constrCoeff=1;
sigmaNoisePoses=0.05;
sigmaNoiseNormals=0.05;
sigmaNoiseDistances=0.05;

NPlanes=20;

allDatasets=cell(NTrials,1);
allResults=cell(NTrials,1);
allErrorsAvg=cell(NTrials,1);
allErrorsAvgMeasurement=cell(NTrials,1);
allErrorsLinear=cell(NTrials,1);
for iTrial=1:NTrials
    fprintf('Trial %d/%d: ',iTrial,NTrials)
    
    %generate dataset
    [dataset.GCamera,dataset.GDepth,dataset.GPlaneTags,dataset.NVecPlanesDepth]=...
        cameraDepthSensorPlaneCalibration_buildDataset(NPlanes,GTag);

    %add noise
    dataset.GPlaneTagsNoise=noiserigid(dataset.GPlaneTags,sigmaNoisePoses);
    dataset.NVecPlanesDepthNoise=planeRandn(dataset.NVecPlanesDepth,sigmaNoiseNormals,sigmaNoiseDistances);

    %average
    [result.GCameraAvg,result.GDepthAvg,result.GPlaneTagsAvg,result.GDepthLinear]=...
        cameraDepthSensorPlaneCalibration(dataset.GPlaneTagsNoise,dataset.NVecPlanesDepthNoise,...
        'coeffPlaneConstraints',constrCoeff,...
        'optsLocalization',{'optsLieMinimize',{'maxIt',NIt,'progressBar'}});

    %compute errors
    errorAvg=computeError(result.GDepthAvg,result.GPlaneTagsAvg,dataset.NVecPlanesDepthNoise,dataset.GDepth);
    errorAvgMeasurement=computeError(result.GDepthAvg,dataset.GPlaneTagsNoise,dataset.NVecPlanesDepthNoise,dataset.GDepth);
    errorLinear=computeError(result.GDepthLinear,dataset.GPlaneTagsNoise,dataset.NVecPlanesDepthNoise,dataset.GDepth);
    
    %store results of the trial
    allDatasets{iTrial}=dataset;
    allResults{iTrial}=result;
    allErrorsAvg{iTrial}=errorAvg;
    allErrorsAvgMeasurement{iTrial}=errorAvgMeasurement;
    allErrorsLinear{iTrial}=errorLinear;
end

save([mfilename '_Trials' num2str(NTrials) '_It' num2str(NIt)]);

function error=computeError(GDepthMeasured,GPlaneTagsMeasured,NVecPlanesDepth,GDepthTruth)
[error.residualRot,error.residualTransl]=cameraDepthSensorPlaneCalibration_residualsPlanes(GDepthMeasured,GPlaneTagsMeasured,NVecPlanesDepth);
error.poseRot=rot_dist(G2R(GDepthMeasured),G2R(GDepthTruth));
error.poseTransl=norm(G2T(GDepthMeasured)-G2T(GDepthTruth));
