function localization_MLE_rigid_test
fileName=[mfilename '_' regexprep(datestr(now),'[ :]','_')];
diary([fileName '.txt'])
NTrials=50;
results=repmat(struct('t_node',[],'rotErr',[],'translErr',[],'errors',[],'t_nodeMod',[],'rotErrMod',[],'translErrMod',[],'errorsMod',[]),1,NTrials);
for iTrial=1:NTrials
    t_node=buildTestTagNetwork('seed',[]);
    [t_node,errors]=localization_MLE_rigid(t_node);
    disp('Errors with covariances')
    t_node=testNetworkCompensate(t_node);
    [rotErr,translErr]=testNetworkComputeErrors(t_node,'references');
    rotErr=rotErr*180/pi;
    translErr=translErr*180/pi;
    disp(['Mean rot error    ' num2str(mean(rotErr))])
    disp(['Std  rot error    ' num2str(std(rotErr))])
    disp(['Mean transl error ' num2str(mean(translErr))])
    disp(['Std  transl error ' num2str(std(translErr))])
    %testNetworkShowErrors(t_node,'references')
    results(iTrial).t_node=t_node;
    results(iTrial).rotErr=rotErr;
    results(iTrial).translErr=translErr;
    results(iTrial).errors=errors;

    t_nodeMod=testNetworkAddDispersionMatricesRT(t_node,'methodR','identity','methodT','identity','methodCoupling','zero');
    [t_nodeMod,errorsMod]=localization_MLE_rigid(t_nodeMod);
    disp('Errors without covariances')
    t_nodeMod=testNetworkCompensate(t_nodeMod);
    [rotErrMod,translErrMod]=testNetworkComputeErrors(t_nodeMod,'references');
    rotErrMod=rotErrMod*180/pi;
    translErrMod=translErrMod*180/pi;
    disp(['Mean rot error    ' num2str(mean(rotErrMod))])
    disp(['Std  rot error    ' num2str(std(rotErrMod))])
    disp(['Mean transl error ' num2str(mean(translErrMod))])
    disp(['Std  transl error ' num2str(std(translErrMod))])
    %testNetworkShowErrors(t_nodeMod,'references')
    results(iTrial).t_nodeMod=t_nodeMod;
    results(iTrial).rotErrMod=rotErrMod;
    results(iTrial).translErrMod=translErrMod;
    results(iTrial).errorsMod=errorsMod;
    save(fileName)
end
diary off
