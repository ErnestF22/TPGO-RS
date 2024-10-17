function multiCalibrationNew_testRealData
resetRands();
NTrials=30;

disp('# Loading correspondences')
%correspondences=multiCalibration_datasetToCorrespondencesFile('calibrationDataset3SensorsSubset');
correspondences=multiCalibration_datasetToCorrespondencesFile('calibrationDataset7Sensors');
correspondences=multiCalibration_correspondencesProcess(correspondences);

errors=cell(NTrials,1);
%test with multiple random choices for the initialization tree
for iTrial=1:NTrials
    disp(['# Trial ' num2str(iTrial) '/' num2str(NTrials)])
    
    [posesAbsolute,posesRelativeMeasured,posesAbsoluteTree]=multiCalibrationNew(correspondences);
    posesRelativeAvg=multiCalibration_correspondencesPoses(correspondences,posesAbsolute);
    posesRelativeTree=multiCalibration_correspondencesPoses(correspondences,posesAbsoluteTree);
    errors{iTrial}.averaged=multiCalibration_errorsCompute(posesRelativeAvg,correspondences);
    errors{iTrial}.measured=multiCalibration_errorsCompute(posesRelativeMeasured,correspondences);
    errors{iTrial}.tree=multiCalibration_errorsCompute(posesRelativeTree,correspondences);
end

figure(1)
multiCalibration_posesDisplay(posesAbsolute)

figure(2)
disp('# Errors')
errorsAggregated=multiCalibration_errorsAggregate(errors);
multiCalibration_errorsDisplay(errorsAggregated,correspondences)

save([mfilename '_data'])
