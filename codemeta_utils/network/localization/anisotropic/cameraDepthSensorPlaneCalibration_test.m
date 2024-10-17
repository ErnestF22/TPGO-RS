function cameraDepthSensorPlaneCalibration_test
resetRands(1)
flagNoisy=true;

L=3;
offset=[L/2;-L/2;0];
GTag1=RT2G(eye(3),zeros(3,1)+offset);
GTag2=RT2G(eye(3),[-L;0;0]+offset);
GTag3=RT2G(eye(3),[0;L;0]+offset);
GTag=cat(3,GTag1,GTag2,GTag3);
NIt=3000;
constrCoeff=1;

NPlanes=20;

methodAbsolutePoses='references';

if ~flagNoisy
    sigmaNoisePoses=0;
    sigmaNoiseNormals=0;
    sigmaNoiseDistances=0;
else
    sigmaNoisePoses=0.05;
    sigmaNoiseNormals=0.05;
    sigmaNoiseDistances=0.05;
end

%generate dataset
[GCamera,GDepth,GPlaneTags,NVecPlanesDepth]=...
    cameraDepthSensorPlaneCalibration_buildDataset(NPlanes,GTag);

%add noise
GPlaneTagsNoise=noiserigid(GPlaneTags,sigmaNoisePoses);
NVecPlanesDepthNoise=planeRandn(NVecPlanesDepth,sigmaNoiseNormals,sigmaNoiseDistances);

%average
[GCameraAvg,GDepthAvg,GPlaneTagsAvg,GDepthLinear]=...
    cameraDepthSensorPlaneCalibration(GPlaneTagsNoise,NVecPlanesDepthNoise,...
    'coeffPlaneConstraints',constrCoeff,...
    'optsLocalization',{'optsLieMinimize',{'maxIt',NIt,'progressBar'}});

save([mfilename '_data_It' num2str(NIt) '_coeff' num2str(constrCoeff) '_Planes' num2str(NPlanes) '.mat'])

%show results
NVecPlanes=rigidTransformG(GDepth,NVecPlanesDepth,'planes',methodAbsolutePoses,'cw');

draw3dcameraFromG(GPlaneTagsAvg,methodAbsolutePoses,'shape','april')
hold on
draw3dcameraFromG(GCameraAvg,methodAbsolutePoses,'color','g')
draw3dcameraFromG(GDepthAvg,methodAbsolutePoses,'color','b')
draw3dPlane(NVecPlanes,'side',3)
axis equal

disp('Error GDepthAvg')
displayPoseError(GDepthAvg,GDepth)
[eRot,eTransl]=cameraDepthSensorPlaneCalibration_residualsPlanes(GDepthAvg,GPlaneTagsAvg,NVecPlanesDepthNoise);
disp(['Mean normal angle residual [deg]: ' num2str(mean(eRot)*180/pi)])
disp(['Mean absolute translation residual [cm]: ' num2str(mean(abs(eTransl))*100)])
disp('Error GDepthLinear')
displayPoseError(GDepthLinear,GDepth)
[eRot,eTransl]=cameraDepthSensorPlaneCalibration_residualsPlanes(GDepthLinear,GPlaneTagsNoise,NVecPlanesDepthNoise);
disp(['Mean normal angle residual [deg]: ' num2str(mean(eRot)*180/pi)])
disp(['Mean absolute translation residual [cm]: ' num2str(mean(abs(eTransl))*100)])

function displayPoseError(GDepthAvg,GDepth)
disp('\tRotation (relative angle) [deg]')
disp(rot_dist(G2R(GDepthAvg),G2R(GDepth))*180/pi)
disp('\tTranslation (norm of difference) [cm]')
disp(norm(G2T(GDepthAvg)-G2T(GDepth))*100)
