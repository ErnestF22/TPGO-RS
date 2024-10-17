%function [cTrain,cTest]=multiCalibration_correspondencesRawDataSplit(c,fractionTest)
%Split the rawData field in each element of the correspondences cell array
%into training and test sets according to the fractionTest ratio.
%All the other fields are copied verbatim. Before using the correspondences
%in multiCalibration.m, you should call multiCalibration_correspondencesProcess
function [cTrain,cTest]=multiCalibration_correspondencesRawDataSplit(c,fractionTest)

cTrain=c;
cTest=c;

NCorrespondences=length(c);
for iCorrespondence=1:NCorrespondences
    cCurrent=c{iCorrespondence};
    NData=size(cCurrent.rawData{1},1);
    NDataTest=round(fractionTest*NData);
    idxTest=randperm(NData,NDataTest);
    idxTrain=setdiff(1:NData,idxTest);
    
    cTrain{iCorrespondence}.rawData=cellfun(@(x) x(idxTrain,:),cCurrent.rawData,'UniformOutput',false);
    cTest{iCorrespondence}.rawData=cellfun(@(x) x(idxTest,:),cCurrent.rawData,'UniformOutput',false);
end
