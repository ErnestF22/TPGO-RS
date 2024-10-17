function POCMultiCalibration_correspSplit
load multiCalibration_testRealData_data
fractionTest=0.2;

c=correspondences;
fieldNames=fields(c);
for iField=1:length(fieldNames)
    fn=fieldNames{iField};
    [cTrain.(fn),cTest.(fn)]=structSplit(c.(fn),fractionTest);
end

keyboard


function [sTrain,sTest]=structSplit(s,fractionTest)
fieldNames=fields(s);
NData=size(s.(fieldNames{1}),2);
NDataTest=round(fractionTest*NData);
idxTest=randperm(NData,NDataTest);
idxTrain=setdiff(1:NData,idxTest);

for iField=1:length(fieldNames)
    fn=fieldNames{iField};
    if size(s.(fn),2)==NData
        sTrain.(fn)=s.(fn)(:,idxTrain);
        sTest.(fn)=s.(fn)(:,idxTest);
    end
end
