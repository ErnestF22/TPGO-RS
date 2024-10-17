%Extract errors stored in a cell array of structs
function [methodNames,data]=sfm_utilityExtractErrors(errors)

[NTrials,NSigmas,NDatasets]=size(errors);
methodNames=fields(errors{1,1});
NMethods=length(methodNames);

disp('# Statistics:')
fprintf('\t%d datasets\n',NDatasets)
fprintf('\t%d methods\n',NMethods)
fprintf('\t%d sigmas\n',NSigmas)
fprintf('\t%d trials\n',NTrials)

data=cell(1,NDatasets);
for iDataset=1:NDatasets
    dimData=length(errors{1,1,iDataset}.(methodNames{1}));
    fprintf('- Dataset %d: %d edges\n',iDataset,dimData)
    data{iDataset}=zeros(NTrials,NSigmas,dimData,NMethods);
    for iMethod=1:NMethods
        for iTrial=1:NTrials
            for iSigma=1:NSigmas
                data{iDataset}(iTrial,iSigma,:,iMethod)=...
                    errors{iTrial,iSigma,iDataset}.(methodNames{iMethod});
            end
        end
    end
    data{iDataset}=permute(data{iDataset},[2 4 3 1]);
end
