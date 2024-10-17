function multiCalibration_compareErrors
s=load('multiCalibration_testRealData_N30_data');
sCov=load('multiCalibration_testRealData_N30_covariances_data');

figure(1)
subplot(2,2,1)
cumDistBoxPerc([s.errorsAggregated.cameraVelodyne.rot sCov.errorsAggregated.cameraVelodyne.rot]*180/pi)
legend('Isotropic','Anisotropic')
title('Camera-Velodyne rotation')

subplot(2,2,2)
cumDistBoxPerc(abs([s.errorsAggregated.cameraVelodyne.transl sCov.errorsAggregated.cameraVelodyne.transl])*100)
title('Camera-Velodyne translation')

subplot(2,2,3)
cumDistBoxPerc(abs([s.errorsAggregated.hokuyoCamera sCov.errorsAggregated.hokuyoCamera])*100)
title('Hokuyo-Camera')

subplot(2,2,4)
cumDistBoxPerc(abs([s.errorsAggregated.hokuyoVelodyne sCov.errorsAggregated.hokuyoVelodyne])*100)
title('Hokuyo-Velodyne')

