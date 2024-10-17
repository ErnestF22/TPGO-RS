function cvpr17_plotPrecisionRecall(precision,recall,iThreshold)
methodNames=fields(precision);
NMethods=length(methodNames);
for iMethod=1:NMethods
    method=methodNames{iMethod};
    plot(precision.(method)(:,iThreshold),recall.(method)(:,iThreshold),'.-',...
        'linewidth',1,'MarkerSize',10)
    hold on
end
hold off
legend(methodNames)
xlabel('count/nMatches')
ylabel('count/nFeatures')
