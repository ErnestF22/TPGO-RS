function bearingNetworkPreconditioner_test_trials_plot
%load bearingNetworkPreconditioner_test_trials_09-Oct-2014_09_55_40_processed
load bearingNetworkPreconditioner_test_trials_10-Oct-2014_10_35_25_processed

idxBaseline=[1];
idxAll=1:13;
idxNeumannOpt=[10 11 12 13];
idxNeumannStd=[6 7 8 9];
idxNeumann=[idxNeumannStd idxNeumannOpt];
idxNeumannLow=[6 7 10 11];
idxNeumannHigh=[8 9 12 13];
idxSPD=[3 4];
idxBest=[3 11 12];
idx=[idxBest idxBaseline];

NFieldsNames=fieldnames(dataProcessed);
NFieldsNames=NFieldsNames([5 1 4 2 3]);
NNFieldsNames=length(NFieldsNames);
NCols=3;
NRows=ceil(NNFieldsNames/NCols);
figure(1)
for iField=1:NNFieldsNames
    fieldName=NFieldsNames{iField};
    
    subplot(NRows,NCols,iField)
    semilogy(dataProcessed.(fieldName).t,mean(dataProcessed.(fieldName).phi(:,idx,:),3));
    legend(dataProcessed.(fieldName).types(idx))
    title(fieldName)
end
