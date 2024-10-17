function sfm_poseCombine_testTrialsOutliers_plots(idxSelect)
sigmaOutliers=0:0.05:1;
load('sfm_poseCombine_testTrialsOutliers_17-Feb-2015_16_02_52')
%load('sfm_poseCombine_testTrialsOutliers_02-Feb-2015_15_26_39.mat')
%load('sfm_poseCombine_testTrialsOutliers_16-Jan-2015_12_44_04.mat')

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
    plot(sigmaOutliers,mean(e(:,idxSelect,:),3))
    plot(sigmaOutliers,median(e(:,idxSelect,:),3),'--')
    hold off
    legend(methodNames(idxSelect))
    xlabel('Fraction of outliers')
    ylabel('Edge average error [deg]')
end

