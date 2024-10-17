function bearingNetworkPreconditioner_test_trials
fileNameSave=[mfilename '_' regexprep(datestr(now),'[ :]','_')];

NNodes=[7 12 15 10 5];
TDelay=[5 13 15 10 5];
NTrials=100;

data=struct();

for iNodes=1:length(NNodes)
    disp(['# Nodes ' num2str(iNodes) '/' num2str(NNodes)])
    outputTrials=cell(1,NTrials);
    for iTrial=1:NTrials
        disp(['# Trial ' num2str(iTrial) '/' num2str(NTrials)])
        outputTrials{iTrial}=bearingNetworkPreconditioner_test_singleTrial(NNodes(iNodes),TDelay(iNodes));
        if mod(iNodes,15)==1
            data.(['N' num2str(NNodes(iNodes))])=outputTrials(1:iTrial);
            save(fileNameSave,'data')
        end
    end
    data.(['N' num2str(NNodes(iNodes))])=outputTrials;
    save(fileNameSave,'data')
end
        