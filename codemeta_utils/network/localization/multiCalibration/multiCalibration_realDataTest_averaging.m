function multiCalibration_realDataTest_averaging
load multiCalibration_realDataTest_data
[t_node,idxCamera,idxVelodyne,idxHokuyo]=multiCalibrationPrepareNetwork(GCameraVelodyneRansac,GHokuyoCameraRansac,GHokuyoVelodyneRansac);
optsLocalization={'displayIt','showMessages'};
t_node=localization_MLE_rigid(t_node,optsLocalization{:},'noinit');
% Compute poses
methodAbsolutePoses='reference';
figure(2)
testNetworkDisplayErrors(t_node,'rt','methodAbsolutePoses',methodAbsolutePoses,'methodTranslationError','norm')
figure(1)
draw3dcameraFromG(t_node.gi(:,:,idxCamera),'scale',0.25,'methodAbsolutePoses',methodAbsolutePoses);
hold on
draw3dcameraFromG(t_node.gi(:,:,idxVelodyne),'shape','velodyne','scale',0.25,'methodAbsolutePoses',methodAbsolutePoses);
draw3dcameraFromG(t_node.gi(:,:,idxHokuyo),'shape','hokuyo','scale',0.25,'methodAbsolutePoses',methodAbsolutePoses);
hold off
axis equal

GCamera=t_node.gi(:,:,idxCamera);
GVelodyne=t_node.gi(:,:,idxVelodyne);
GHokuyo=t_node.gi(:,:,idxHokuyo);

save multiCalibration_averaging_result GCamera GVelodyne GHokuyo
