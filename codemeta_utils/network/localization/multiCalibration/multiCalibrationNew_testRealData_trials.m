function multiCalibrationNew_testRealData_trials
resetRands()
flagUseCovariances=false;
NTrials=100;

fileNameData=[mfilename '_N' num2str(NTrials)];
if flagUseCovariances
    fileNameData=[fileNameData '_covariances'];
end

disp('# Loading correspondences')
correspondences=multiCalbration_datasetToCorrespondencesFile('calibrationDataset');
correspondences=multiCalibration_correspondencesProcess(correspondences);


errors=cell(NTrials,1);
poses=cell(NTrials,1);
times=zeros(NTrials,2);
for iTrial=1:NTrials
    disp(['## Trial ' num2str(iTrial) '/' num2str(NTrials)])
    [correspondencesTrain,correspondencesTest]=multiCalibration_correspondencesRawDataSplit(...
        correspondences,0.2);
    
    correspondencesTrain=multiCalibration_correspondencesProcess(correspondencesTrain);
    correspondencesTest=multiCalibration_correspondencesProcess(correspondencesTest);
    
    times(iTrial,1)=cputime;
    [posesAbsoluteAvg,posesRelativeMeasured,posesAbsoluteTree]=multiCalibrationNew(correspondencesTrain,...
        'flagUseCovariances',flagUseCovariances);
    times(iTrial,2)=cputime;

    posesRelativeAvg=multiCalibration_correspondencesPoses(correspondencesTrain,posesAbsoluteAvg);
    posesRelativeTree=multiCalibration_correspondencesPoses(correspondences,posesAbsoluteTree);
    errors{iTrial}.measured=multiCalibration_errorsCompute(posesRelativeMeasured,correspondencesTest);
    errors{iTrial}.averaged=multiCalibration_errorsCompute(posesRelativeAvg,correspondencesTest);
    errors{iTrial}.tree=multiCalibration_errorsCompute(posesRelativeTree,correspondencesTest);
    
    %save([fileNameData '_data_partial'])
end

save([fileNameData '_data'])

figure(1)
multiCalibration_posesDisplay(posesAbsoluteAvg)

figure(2)
disp('# Errors')
errorsAggregated=multiCalibration_errorsAggregate(errors);
multiCalibration_errorsDisplay(errorsAggregated,correspondences)


