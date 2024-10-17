function sfm_poseCombine_testTrials_plots(idxSelect)
sigmaNoise=0.001*(0:2.5:20);
load('sfm_poseCombine_testTrials_17-Feb-2015_15_58_08')
%load('sfm_poseCombine_testTrials_16-Jan-2015_12_44_09')


[methodNames,data]=sfm_utilityExtractErrors(errors);
if ~exist('idxSelect','var')
    idxSelect=1:length(methodNames);
end
colors=rbg(length(idxSelect));

NDatasets=length(data);
for iDataset=1:NDatasets
    figure(iDataset)
    sz=size(data{iDataset});
    e=reshape(data{iDataset},sz(1),sz(2),[])*180/pi;
    cla
    set(gca,'ColorOrder',colors);
    hold on
    plot(sigmaNoise,mean(e(:,idxSelect,:),3))
    plot(sigmaNoise,median(e(:,idxSelect,:),3),'--')
    hold off
    legend(methodNames(idxSelect))
    xlabel('Noise level (variance in normalized coords)')
    ylabel('Edge average error [deg]')
end
