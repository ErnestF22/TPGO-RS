function c=multiCalibration_preprocessCorrespondences(c)
hc=c.hokuyoCamera;
hv=c.hokuyoVelodyne;

[hv.normalVelodyneFlat,hv.distanceVelodyneFlat,hv.pointsHokuyo2D]=prepareDataHokuyo(hv.normalVelodyne,hv.distanceVelodyne,hv.pointsHokuyo);
[hc.normalCameraFlat,hc.distanceCameraFlat,hc.pointsHokuyo2D]=prepareDataHokuyo(hc.normalCamera,hc.distanceCamera,hc.pointsHokuyo);

c.hokuyoCamera=hc;
c.hokuyoVelodyne=hv;

function [normalCameraFlat,distanceCameraFlat,pointsHokuyo2D]=prepareDataHokuyo(normalCamera,distanceCamera,pointsHokuyo)
normalCameraFlat=reshape(repmat(normalCamera,[2 1]),3,[]);
distanceCameraFlat=reshape(repmat(distanceCamera,[2 1]),1,[]);
pointsHokuyo2D=reshape(pointsHokuyo([2 3 5 6],:),2,[]);
