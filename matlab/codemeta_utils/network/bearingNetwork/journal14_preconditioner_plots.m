function journal14_preconditioner_plots
%generates plots with cost for preconditioner in the journal paper
figDir='../../../papers/network/TRO14-BearingFormation/figures/';
load bearingNetworkPreconditioner_test_trials_10-Oct-2014_10_35_25_processed
flagFigSave=2;
figDim=[250 135];

fs=setFigFontSize(8);
fn=setFigFont('Times');


idxBaseline=[1];
idxBest=[3 11 12];
idx=[idxBaseline idxBest];

NFieldsNames=fieldnames(dataProcessed);
NFieldsNames=NFieldsNames([5 1 4 2 3]);
NNFieldsNames=length(NFieldsNames);
switchTimes=[5 5 10 13 15];
figure(1)
for iField=1:NNFieldsNames
    fieldName=NFieldsNames{iField};
    
    figFile=['preconditioner_cost_' fieldName];
    
    semilogy(dataProcessed.(fieldName).t,mean(dataProcessed.(fieldName).phi(:,idx,:),3));
    if iField==NNFieldsNames-1
        %legend(dataProcessed.(fieldName).types(idx))
        legend('No acceleration','Min. Frobenius norm','Neumann series, p=2', 'Neumann series, p=10')
    end
    
    ax=axis;
    hold on
    plot(switchTimes(iField)*[1 1],[ax(3) ax(4)],'k--');
    hold off
    
    savefigure(fullfile(figDir,figFile),'epsc',figDim,flagFigSave)
end

setFigFontSize(fs);
setFigFont(fn);

