function localization_MLE_rigid_testTrials
flagUseIsotropicInit=true;
fileName=[mfilename '_' regexprep(datestr(now),'[ :]','_') '_I' num2str(flagUseIsotropicInit)];
diary([fileName '.txt'])
NTrials=50;
results=repmat(struct('isotropic',[],'weighted',[],'covariances',[],'spectral',[]),NTrials,1);
for iTrial=1:NTrials
    disp(['Trial ' num2str(iTrial) '/' num2str(NTrials)])
    t_node=buildTestTagNetwork('seed',[]);
    result(iTrial)=localization_MLE_rigid_testSingle(t_node,'flagUseIsotropicInit',flagUseIsotropicInit);
    save(fileName)
end
diary off
