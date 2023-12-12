function bearingNetworkPreconditioner_test_trials_process
%fileName='bearingNetworkPreconditioner_test_trials_09-Oct-2014_09_55_40';
fileName='bearingNetworkPreconditioner_test_trials_10-Oct-2014_10_35_25';
load([fileName '.mat'])
dataProcessed=struct();
NFieldsNames=fieldnames(data);
NNFields=length(NFieldsNames);

for iNFields=1:NNFields
    NFieldName=NFieldsNames{iNFields};
    typeFields=fieldnames(data.(NFieldName){1}.t);
    NTypeFields=length(typeFields);
    
    NTrials=length(data.(NFieldName));
    tMax=max(data.(NFieldName){1}.t.(typeFields{1}));
    t=0:0.1:tMax;
    Nt=length(t);
    
    phi=zeros(Nt,NTypeFields,NTrials);
    for iTrial=1:NTrials
        for iType=1:NTypeFields
            type=typeFields{iType};
            t1=data.(NFieldName){iTrial}.t.(type);
            phi1=data.(NFieldName){iTrial}.phi.(type);
            phi(:,iType,iTrial)=interp1(t1,phi1,t);
        end
    end
    
    dataProcessed.(NFieldName).t=t;
    dataProcessed.(NFieldName).phi=phi;
    dataProcessed.(NFieldName).types=typeFields;
end
 
save([fileName '_processed'],'dataProcessed')

