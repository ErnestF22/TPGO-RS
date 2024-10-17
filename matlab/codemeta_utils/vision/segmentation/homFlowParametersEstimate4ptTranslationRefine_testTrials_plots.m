function homFlowParametersEstimate4ptTranslationRefine_testTrials_plots
load('homFlowParametersEstimate4ptRefine_testTrials_data')

errorsMean=zeros(4,NSigmas);
legendText=cell(2,1);
for iSigma=1:NSigmas
    errorsMean(1,iSigma)=mean([errors(:,iSigma).linear]);
    legendText{1}='Linear';
    errorsMean(2,iSigma)=mean([errors(:,iSigma).refined2]);
    legendText{2}='Refined 2';
    errorsMean(3,iSigma)=mean([errors(:,iSigma).refined20]);
    legendText{3}='Refined 20';
    errorsMean(4,iSigma)=mean([errors(:,iSigma).refined50]);
    legendText{4}='Refined 50';
end

plot(sigma,errorsMean)
legend(legendText)

