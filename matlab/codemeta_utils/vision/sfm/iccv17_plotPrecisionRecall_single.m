function iccv17_plotPrecisionRecall_single(precision,recall,iThreshold)
lineOpts={'.-','linewidth',1,'MarkerSize',10};
plot(precision.Pair(:,iThreshold),recall.Pair(:,iThreshold),lineOpts{:},'Color','k')
hold on
plot(precision.QuickMatch(:,iThreshold),recall.QuickMatch(:,iThreshold),lineOpts{:},'Color','b')
hold off
xlim([0 1]); set(gca,'xtick',0:0.2:1);
grid on
legend('Pairwise','QuickMatch')
xlabel('count/nMatches')
ylabel('count/nFeatures')
