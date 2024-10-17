function errors=multiCalibration_computeErrors(poses,correspondences)
c=correspondences;
cv=c.cameraVelodyne;
hc=c.hokuyoCamera;
hv=c.hokuyoVelodyne;

%define shorthands for functions for computing residuals
errCV=@(G) poseEstimationFromNormalsDistancesResidualsG(...
    G,...
    cv.normalVelodyne,cv.distanceVelodyne,...
    cv.normalCamera,cv.distanceCamera);
errHC=@(G) poseEstimationFromNormalsPointsResidualsG(...
    G,hc.normalCameraFlat,hc.distanceCameraFlat,hc.pointsHokuyo2D);
errHV=@(G) poseEstimationFromNormalsPointsResidualsG(...
    G,hv.normalVelodyneFlat,hv.distanceVelodyneFlat,hv.pointsHokuyo2D);

%collect residuals
GCameraVelodyne=computeRelativePoseFromG(poses.GVelodyne,poses.GCamera,'references');
[errors.cameraVelodyneRot,errors.cameraVelodyneTransl]=...
    errCV(GCameraVelodyne);
[errors.cameraVelodyneRotMeasured,errors.cameraVelodyneTranslMeasured]=...
    errCV(poses.GCameraVelodyne);
[errors.cameraVelodyneRotHand,errors.cameraVelodyneTranslHand]=...
    errCV(poses.GCameraVelodyneHand);
GCameraVelodyneIndirect=poses.GHokuyoVelodyne*invg(poses.GHokuyoCamera);
[errors.cameraVelodyneRotIndirect,errors.cameraVelodyneTranslIndirect]=...
    errCV(GCameraVelodyneIndirect);


GHokuyoCamera=computeRelativePoseFromG(poses.GCamera,poses.GHokuyo,'references');
GHokuyoCameraIndirect=invg(poses.GCameraVelodyne)*poses.GHokuyoVelodyne;
errors.hokuyoCamera=errHC(GHokuyoCamera);
errors.hokuyoCameraMeasured=errHC(poses.GHokuyoCamera);
errors.hokuyoCameraHand=errHC(poses.GHokuyoCameraHand);
errors.hokuyoCameraIndirect=errHC(GHokuyoCameraIndirect);

GHokuyoVelodyne=computeRelativePoseFromG(poses.GVelodyne,poses.GHokuyo,'references');
GHokuyoVelodyneIndirect=poses.GCameraVelodyne*poses.GHokuyoCamera;
errors.hokuyoVelodyne=errHV(GHokuyoVelodyne);
errors.hokuyoVelodyneMeasured=errHV(poses.GHokuyoVelodyne);
errors.hokuyoVelodyneHand=errHV(poses.GHokuyoVelodyneHand);
errors.hokuyoVelodyneIndirect=errHV(GHokuyoVelodyneIndirect);
