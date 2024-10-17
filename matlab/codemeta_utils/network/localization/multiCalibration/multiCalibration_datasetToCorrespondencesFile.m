function c=multiCalibration_datasetToCorrespondencesFile(dirName)
%the most important fields are:
%nodeNames (used for building edges)
%rawData (contains the data from the files)
%type (used for knowing how to use the rawData)

fileName=fullfile(dirName,'full_cal_info.mat');
if ~exist(fileName,'file')
    error('Cannot access full_cal_info.mat to retrive node names')
end
s=load(fileName);
nodeNames=fieldnames(s);
d=dir(fullfile(dirName,'*.mat'));
NFiles=length(d);
c={};
for iFiles=1:NFiles
    fileName=fullfile(dirName,d(iFiles).name);
    
    [idxPattern,idxChar]=multiStrFind(fileName,nodeNames);
    
    if isempty(idxPattern)
        disp(['Skipped ' fileName])
        continue
    end
    if length(idxPattern)~=2
        warning([fileName ' has an ambiguous name interpretation'])
        continue
    end
    
    nodeNamesFile=nodeNames(idxPattern);

    %load raw data
    s=load(fileName);
    rawData=cell(1,2);
    for iData=1:2
        fieldName=[nodeNamesFile{iData} '_detections'];
        if ~isfield(s,fieldName)
            error([fileName ' does not contain expected fields'])
        end
        rawData{iData}=s.(fieldName);
    end
    
    %detect if the file contains points or only planes
    flagHasPoints=cellfun(@(x) ~isempty(x),strfind(nodeNamesFile,'hokuyo'));
    if any(flagHasPoints)
        type='plane-point';
        %make sure that point data are in the second element of rawData
        if flagHasPoints(1)
            rawData=fliplr(rawData);
            nodeNamesFile=flipud(nodeNamesFile);
        end
    else
        type='plane-plane';
    end
    
    cCurrent=[];
    cCurrent.fileName=fileName;
    cCurrent.nodeNames=nodeNamesFile;
    cCurrent.type=type;
    cCurrent.rawData=rawData;
    
    c{end+1}=cCurrent;
end

%Runs strfind(str,patterns{i}) for all i and records which patterns have
%been found and at what character
function [idxPattern,idxChar]=multiStrFind(str,patterns)
idxPattern=[];
idxChar=[];
for idxP=1:length(patterns)
    idxC=strfind(str,patterns{idxP});
    if ~isempty(idxC)
        idxPattern=[idxPattern idxP];
        idxChar=[idxChar idxC];
    end
end