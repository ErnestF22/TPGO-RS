function multiCalibration_realData_process
load calibrationDataset/april_vel.mat
[normalCamera,distanceCamera]=cnormalize(apriltag_planes');
distanceCamera=-distanceCamera;
[normalVelodyne,distanceVelodyne]=cnormalize(velodyne_planes');
distanceVelodyne=-distanceVelodyne;
save calibrationDataset/april_vel_processed.mat normalCamera normalVelodyne distanceCamera distanceVelodyne

load calibrationDataset/april_hok.mat
[normalCamera,distanceCamera,normalCameraFlat,distanceCameraFlat,hokuyo_endpts2D]=prepareDataCameraHokuyo(apriltag_planes,hokuyo_endpts);
save calibrationDataset/april_hok_processed.mat normalCamera distanceCamera normalCameraFlat distanceCameraFlat hokuyo_endpts2D

load calibrationDataset/hok_vel.mat
[normalVelodyne,distanceVelodyne,normalVelodyneFlat,distanceVelodyneFlat,hokuyo_endpts2D]=prepareDataCameraHokuyo(velodyne_planes,hokuyo_endpts);
save calibrationDataset/hok_vel_processed.mat normalVelodyne distanceVelodyne normalVelodyneFlat distanceVelodyneFlat hokuyo_endpts2D


function [normalCamera,distanceCamera,normalCameraFlat,distanceCameraFlat,hokuyo_endpts2D]=prepareDataCameraHokuyo(apriltag_planes,hokuyo_endpts)
[normalCamera,distanceCamera]=cnormalize(apriltag_planes');
distanceCamera=-distanceCamera;
normalCameraFlat=reshape(repmat(normalCamera,[2 1]),3,[]);
distanceCameraFlat=reshape(repmat(distanceCamera,[2 1]),1,[]);
hokuyo_endpts2D=reshape(hokuyo_endpts(:,[2 3 5 6])',2,[]);
