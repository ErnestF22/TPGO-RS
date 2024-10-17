function [poses,output]=multiCalibration(correspondences,varargin)
flagUseCovariances=false;
NTrials=1000;
optsRansac={'Ntrials',NTrials,'collectResiduals','inlierEstimate','waitBar'};
c=correspondences;
cv=c.cameraVelodyne;
hc=c.hokuyoCamera;
hv=c.hokuyoVelodyne;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'flagusecovariances'
            ivarargin=ivarargin+1;
            flagUseCovariances=varargin{ivarargin};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

%%% Relative pose estimation 
disp('# Relative pose estimation')

%% Camera-Velodyne calibration
ransacThreshold=0.1;
[poses.GCameraVelodyne,output.outputCameraVelodyne]=poseEstimationFromNormalsDistancesGRansac(...
    cv.normalVelodyne,cv.distanceVelodyne,...
    cv.normalCamera,cv.distanceCamera,ransacThreshold,...
    'optsRansac',optsRansac);
fInliers=output.outputCameraVelodyne.flagInliers;
poses.SigmaCameraVelodyne=poseEstimationFromNormalsDistancesCovarianceG(poses.GCameraVelodyne,...
    cv.normalVelodyne(:,fInliers),cv.distanceVelodyne(fInliers),...
    cv.normalCamera(:,fInliers),cv.distanceCamera(fInliers),...
    0.02,0.001);

%% Camera-Hokuyo calibration
ransacThreshold=0.05;
[poses.GHokuyoCamera,output.outputHokuyoCamera]=poseEstimationFromNormalsPointsGRansac(...
    hc.normalCameraFlat,hc.distanceCameraFlat,hc.pointsHokuyo2D,ransacThreshold,...
    'optsRansac',optsRansac,'refineEnd');

fInliers=output.outputHokuyoCamera.flagInliers;
poses.SigmaHokuyoCamera=poseEstimationFromNormalsPointsCovarianceG(poses.GHokuyoCamera,...
    hc.normalCameraFlat(:,fInliers),hc.distanceCameraFlat(fInliers),hc.pointsHokuyo2D(:,fInliers),...
    0.02,0.02,0.02);

%% Velodyne-Hokuyo calibration
ransacThreshold=0.05;
[poses.GHokuyoVelodyne,output.outputHokuyoVelodyne]=poseEstimationFromNormalsPointsGRansac(...
    hv.normalVelodyneFlat,hv.distanceVelodyneFlat,hv.pointsHokuyo2D,ransacThreshold,...
    'optsRansac',optsRansac,'refineEnd');

fInliers=output.outputHokuyoVelodyne.flagInliers;
poses.SigmaHokuyoVelodyne=poseEstimationFromNormalsPointsCovarianceG(poses.GHokuyoVelodyne,...
    hv.normalVelodyneFlat(:,fInliers),hv.distanceVelodyneFlat(fInliers),hv.pointsHokuyo2D(:,fInliers),...
    0.02,0.02,0.02);

%%% Absolute pose estimation through averaging
disp('# Absolute pose estimation through averaging')
[t_node,idxCamera,idxVelodyne,idxHokuyo]=multiCalibrationPrepareNetwork(poses,flagUseCovariances);
optsLocalization={'showMessages'};
t_node=localization_MLE_rigid(t_node,optsLocalization{:},'noinit');

GCompensate=t_node.gi(:,:,idxCamera);
poses.GCamera=invg(GCompensate)*t_node.gi(:,:,idxCamera);
poses.GVelodyne=invg(GCompensate)*t_node.gi(:,:,idxVelodyne);
poses.GHokuyo=invg(GCompensate)*t_node.gi(:,:,idxHokuyo);
%Note: absolute poses are in 'reference' interpretation
