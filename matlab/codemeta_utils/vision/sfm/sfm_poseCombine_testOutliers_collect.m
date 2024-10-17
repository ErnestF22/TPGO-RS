function sfm_poseCombine_testOutliers_collect(datasets)
if ~exist('datasets','var') || isempty(datasets)
    datasets='castle';
end
if ~iscell(datasets)
    datasets={datasets};
end
dataDir='sfm_poseCombine_testOutliers_data';

eCollected.mean=[];
eCollected.median=[];
eCollected.names={};
eCollected.e={};
eCollected.datasets={};
fileNameSave=fullfile(dataDir,'collected');

for iDataset=1:length(datasets)
    thisDataset=datasets{iDataset};
    fprintf('# Dataset: %s\n',thisDataset);
    d=dir(fullfile(dataDir,[thisDataset '_*']));
    NFiles=length(d);
    s=load(fullfile(dataDir,thisDataset));
    for iFile=1:NFiles
        fileName=d(iFile).name;
        fprintf('Processing: %s\n',fileName);
        s=load(fullfile(dataDir,fileName));
        e=cell2mat(s.errors);
        name=strrep(d(iFile).name,'.mat','');
        name=strrep(name,[thisDataset '_'],'');
        if any(isnan(e(:)))
            warning(['NaN found in data for dataset ' thisDataset ' method ' name])
        end
        [idx,flagNew]=name2idx(name,eCollected.names);
        if flagNew
            eCollected.names{idx}=name;
            eCollected.e{idx}=e;
            eCollected.dataset{idx}={thisDataset};
        else
            eCollected.e{idx}=[eCollected.e{idx}; e];
            eCollected.dataset{idx}=[eCollected.dataset{idx} thisDataset];
        end
    end
    fileNameSave=[fileNameSave '_' thisDataset];
end
sigmaOutliers=s.sigmaOutliers;
NSigmas=length(sigmaOutliers);
NNames=length(eCollected.names);
eCollected.mean=zeros(NSigmas,NNames);
for iName=1:NNames
    e=eCollected.e{iName};
    eCollected.mean(:,iName)=nanmean(e)';
    eCollected.median(:,iName)=nanmedian(eCollected.e{iName})';
end

fprintf('Saving to %s\n',fileNameSave)
save(fileNameSave,'eCollected','sigmaOutliers')

function [idx,flagNew]=name2idx(name,allNames)
flagIdx=strcmp(name,allNames);
if any(flagIdx)
    idx=find(flagIdx,1,'first');
    flagNew=false;
else
    idx=length(allNames)+1;
    flagNew=true;
end
