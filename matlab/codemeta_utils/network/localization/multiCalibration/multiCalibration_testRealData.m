function multiCalibration_testRealData
resetRands()
flagUseCovariances=false;
NTrials=100;

fileNameData=[mfilename '_N' num2str(NTrials)];
if flagUseCovariances
    fileNameData=[fileNameData '_covariances'];
end

load('calibrationDataset/cal_transforms.mat')
load('calibrationDataset/april_vel.mat','apriltag_planes','velodyne_planes')
[normalCamera,distanceCamera]=cnormalize(apriltag_planes');
cv.normalCamera=normalCamera;
cv.distanceCamera=-distanceCamera;
[normalVelodyne,distanceVelodyne]=cnormalize(velodyne_planes');
cv.normalVelodyne=normalVelodyne;
cv.distanceVelodyne=-distanceVelodyne;

load('calibrationDataset/april_hok.mat','apriltag_planes','hokuyo_endpts')
[normalCamera,distanceCamera]=cnormalize(apriltag_planes');
hc.normalCamera=normalCamera;
hc.distanceCamera=-distanceCamera;
hc.pointsHokuyo=hokuyo_endpts';

load('calibrationDataset/hok_vel.mat','velodyne_planes','hokuyo_endpts')
[normalVelodyne,distanceVelodyne]=cnormalize(velodyne_planes');
hv.normalVelodyne=normalVelodyne;
hv.distanceVelodyne=-distanceVelodyne;
hv.pointsHokuyo=hokuyo_endpts';

correspondences.cameraVelodyne=cv;
correspondences.hokuyoCamera=hc;
correspondences.hokuyoVelodyne=hv;

correspondences=multiCalibration_preprocessCorrespondences(correspondences);

errors=cell(NTrials,1);
poses=cell(NTrials,1);
output=cell(NTrials,1);
times=zeros(NTrials,2);
for iTrial=1:NTrials
    disp(['## Trial ' num2str(iTrial) '/' num2str(NTrials)])
    [correspondencesTrain,correspondencesTest]=multiCalibration_splitCorrespondences(...
        correspondences,0.2);

    times(iTrial,1)=cputime;
    [poses{iTrial},output{iTrial}]=multiCalibration(correspondencesTrain,...
        'flagUseCovariances',flagUseCovariances);
    poses{iTrial}.GHokuyoCameraHand=invg(april_to_hok);
    poses{iTrial}.GCameraVelodyneHand=april_to_vel;
    poses{iTrial}.GHokuyoVelodyneHand=hok_to_vel;
    
    errors{iTrial}=multiCalibration_computeErrors(poses{iTrial},correspondencesTest);
    times(iTrial,2)=cputime;
    %save([fileNameData '_data_partial'])
end

methodAbsolutePoses='reference';
figure(1)
draw3dcameraFromG(poses{end}.GCamera,'scale',0.25,'methodAbsolutePoses',methodAbsolutePoses);
hold on
draw3dcameraFromG(poses{end}.GVelodyne,'shape','velodyne','scale',0.25,'methodAbsolutePoses',methodAbsolutePoses);
draw3dcameraFromG(poses{end}.GHokuyo,'shape','hokuyo','scale',0.25,'methodAbsolutePoses',methodAbsolutePoses);
hold off
axis equal

errorsAggregated=multiCalibration_testRealData_aggregateErrors(errors);

disp('CV: Median error normals [deg]')
disp(median(errorsAggregated.cameraVelodyneRot)*180/pi)
disp('CV: Median error distances [cm]')
disp(median(abs(errorsAggregated.cameraVelodyneTransl)*100))
disp('HC: Median error reprojection points [cm]')
disp(median(abs(errorsAggregated.hokuyoCamera)*100))
disp('HV: Median error reprojection points [cm]')
disp(median(abs(errorsAggregated.hokuyoVelodyne)*100))

save([fileNameData '_data'])

