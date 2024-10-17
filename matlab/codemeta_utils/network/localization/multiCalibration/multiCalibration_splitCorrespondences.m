function [correspondencesTrain,correspondencesTest]=multiCalibration_splitCorrespondences(correspondences,fractionTest)

c=correspondences;
fieldNames=fields(c);
for iField=1:length(fieldNames)
    fn=fieldNames{iField};
    [cTrain.(fn),cTest.(fn)]=structSplit(c.(fn),fractionTest);
end

cTrain=multiCalibration_preprocessCorrespondences(cTrain);
cTest=multiCalibration_preprocessCorrespondences(cTest);

correspondencesTrain=cTrain;
correspondencesTest=cTest;

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
