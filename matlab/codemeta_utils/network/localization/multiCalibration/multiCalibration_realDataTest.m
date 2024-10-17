function multiCalibration_realDataTest
resetRands(1)
NTrials=5000;

optsRansac={'Ntrials',NTrials,'collectResiduals','inlierEstimate','waitBar'};
%% Camera-Velodyne calibration
load calibrationDataset/april_vel_processed.mat
GCameraVelodyne=poseEstimationFromNormalsDistancesG(normalVelodyne,distanceVelodyne,normalCamera,distanceCamera);
[GCameraVelodyneRansac,output]=poseEstimationFromNormalsDistancesGRansac(normalVelodyne,distanceVelodyne,normalCamera,distanceCamera,0.1,...
    'optsRansac',optsRansac);
[eRotCameraVelodyne,eTranslCameraVelodyne]=poseEstimationFromNormalsDistancesResidualsG(GCameraVelodyne,normalVelodyne,distanceVelodyne,normalCamera,distanceCamera);
[eRotCameraVelodyneRansac,eTranslCameraVelodyneRansac]=poseEstimationFromNormalsDistancesResidualsG(GCameraVelodyneRansac,normalVelodyne,distanceVelodyne,normalCamera,distanceCamera);
figure(1)
subplot(2,2,1)
cumDistBoxPerc([eRotCameraVelodyne eRotCameraVelodyneRansac]*180/pi)
title('Camera/Velodyne angle residuals')
xlabel('Degrees')
subplot(2,2,2)
cumDistBoxPerc(abs([eTranslCameraVelodyne eTranslCameraVelodyneRansac]))
legend('All data','Ransac')
title('Camera/Velodyne translation residuals')
xlabel('Meters?')
disp('[median(eRotCameraVelodyne) median(eRotCameraVelodyneRansac)]*180/pi');
disp([median(eRotCameraVelodyne) median(eRotCameraVelodyneRansac)]*180/pi);
disp('[median(abs(eTranslCameraVelodyne)) median(abs(eTranslCameraVelodyneRansac))]')
disp([median(abs(eTranslCameraVelodyne)) median(abs(eTranslCameraVelodyneRansac))])

%% Camera-Hokuyo calibration
% load calibrationDataset/april_hok.mat
% [normalCameraFlat,distanceCameraFlat,hokuyo_endpts2D]=prepareDataCameraHokuyo(apriltag_planes,hokuyo_endpts);
load calibrationDataset/april_hok_processed.mat
ransacThreshold=0.05;
GHokuyoCamera=poseEstimationFromNormalsPointsG(normalCameraFlat,distanceCameraFlat,hokuyo_endpts2D);
[GHokuyoCameraRansac,output]=poseEstimationFromNormalsPointsGRansac(normalCameraFlat,distanceCameraFlat,hokuyo_endpts2D,ransacThreshold,...
    'optsRansac',optsRansac);
[GHokuyoCameraRansacRefine,outputRefine]=poseEstimationFromNormalsPointsGRansac(normalCameraFlat,distanceCameraFlat,hokuyo_endpts2D,ransacThreshold,...
    'optsRansac',optsRansac,'refine');
eHokuyoCamera=poseEstimationFromNormalsPointsResidualsG(GHokuyoCamera,normalCameraFlat,distanceCameraFlat,hokuyo_endpts2D);
eHokuyoCameraRansac=poseEstimationFromNormalsPointsResidualsG(GHokuyoCameraRansac,normalCameraFlat,distanceCameraFlat,hokuyo_endpts2D);
eHokuyoCameraRansacRefine=poseEstimationFromNormalsPointsResidualsG(GHokuyoCameraRansacRefine,normalCameraFlat,distanceCameraFlat,hokuyo_endpts2D);
subplot(2,2,3)
cumDistBoxPerc(abs([eHokuyoCamera eHokuyoCameraRansac eHokuyoCameraRansacRefine]))
legend('All data','Ransac')
title('Camera/Hokuyo residuals')
disp('median(abs(eHokuyoCamera))')
disp([median(abs(eHokuyoCamera)) median(abs(eHokuyoCameraRansac)) median(abs(eHokuyoCameraRansacRefine))])

%% Velodyne-Hokuyo calibration
% load calibrationDataset/hok_vel.mat
% [normalVelodyneFlat,distanceVelodyneFlat,hokuyo_endpts2D]=prepareDataCameraHokuyo(velodyne_planes,hokuyo_endpts);
load calibrationDataset/hok_vel_processed.mat
GHokuyoVelodyne=poseEstimationFromNormalsPointsG(normalVelodyneFlat,distanceVelodyneFlat,hokuyo_endpts2D);
[GHokuyoVelodyneRansac,output]=poseEstimationFromNormalsPointsGRansac(normalVelodyneFlat,distanceVelodyneFlat,hokuyo_endpts2D,ransacThreshold,...
    'optsRansac',optsRansac);
[GHokuyoVelodyneRansacRefine,outputRefine]=poseEstimationFromNormalsPointsGRansac(normalVelodyneFlat,distanceVelodyneFlat,hokuyo_endpts2D,ransacThreshold,...
    'optsRansac',optsRansac,'refine');
eHokuyoVelodyne=poseEstimationFromNormalsPointsResidualsG(GHokuyoVelodyne,normalVelodyneFlat,distanceVelodyneFlat,hokuyo_endpts2D);
eHokuyoVelodyneRansac=poseEstimationFromNormalsPointsResidualsG(GHokuyoVelodyneRansac,normalVelodyneFlat,distanceVelodyneFlat,hokuyo_endpts2D);
eHokuyoVelodyneRansacRefine=poseEstimationFromNormalsPointsResidualsG(GHokuyoVelodyneRansacRefine,normalVelodyneFlat,distanceVelodyneFlat,hokuyo_endpts2D);
subplot(2,2,4)
cumDistBoxPerc(abs([eHokuyoVelodyne eHokuyoVelodyneRansac eHokuyoVelodyneRansacRefine]))
legend('All data','Ransac')
title('Velodyne/Hokuyo residuals')
disp('[median(abs(eHokuyoVelodyne)) median(abs(eHokuyoVelodyneRansac))]')
disp([median(abs(eHokuyoVelodyne)) median(abs(eHokuyoVelodyneRansac)) median(abs(eHokuyoVelodyneRansacRefine))])

save([mfilename '_data'])

function [normalCameraFlat,distanceCameraFlat,hokuyo_endpts2D]=prepareDataCameraHokuyo(apriltag_planes,hokuyo_endpts)
[normalCamera,distanceCamera]=cnormalize(apriltag_planes');
normalCameraFlat=reshape(repmat(normalCamera,[2 1]),3,[]);
distanceCameraFlat=reshape(repmat(distanceCamera,[2 1]),1,[]);
hokuyo_endpts2D=reshape(hokuyo_endpts(:,[2 3 5 6])',2,[]);
