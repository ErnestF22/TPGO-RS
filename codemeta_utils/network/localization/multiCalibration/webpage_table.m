function webpage_table
close all
webDataDir='~/Documents/repository/Work/Website/www/assets/';
csvFile=fullfile(webDataDir,'multiCalibrationResults.csv');
figDir=fullfile(webDataDir,'multiCalibrationFigures');
s=load('multiCalibrationNew_testRealData_data');
if exist(csvFile,'file')
    delete(csvFile)
end
diary(csvFile)
multiCalibration_errorsDisplay(s.errorsAggregated,s.correspondences,'save',figDir)
diary off
