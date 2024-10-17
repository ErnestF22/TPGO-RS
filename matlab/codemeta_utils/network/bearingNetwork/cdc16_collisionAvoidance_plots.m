function cdc16_collisionAvoidance_plots
flagPlotTrajectories=false;
%%
close all
figDir='../../../papers/network/cdc16-bearingformation/figures/';
%Column size: 3.25
%Double column size: 7
figDim=[1.3 3/4]*1.7;
figDimSmall=figDim.*[1 0.75];
flagSaveFigure=2;
ax=[-14 8 -8 9.5];
colors=colorDesaturate([1 0 0; 0 1 0; 0 0 1;0 1 1],0.5);
fs=setFigFontSize(8);
fn=setFigFont('Times');

figFileBase=fullfile(figDir,mfilename);
dataFile='cdc16_collisionAvoidance_data';

load(dataFile)

figure(1)
%fake plots to get legend right
for iFlag=1:2
    plot([0,1],[0,1],'-','Visible','off','color',colors(iFlag,:))
    hold on
end
if flagPlotTrajectories
    for iFlag=1:2
        bearingNetworkPlotTrajectories(results{iFlag}.x,'-','color',colors(iFlag,:))
        hold on
    end
end
bearingNetworkPlot(t_node0,'flagPlotCurrent',false)
hold on
for iFlag=1:2
    bearingNetworkPlot(results{iFlag}.t_node,'colorCurrent',colors(iFlag,:),...
        'flagPlotDesired',false);
    hold on
end
hold off
savefigure([figFileBase '_positions'],'pdf',figDim,flagSaveFigure)

figure(2)
for iFlag=1:2
    plot(results{iFlag}.t,min(results{iFlag}.d,[],2),'color',colors(iFlag,:))
    hold on
end
hold off
ax=[0 1500 0 5];
axis(ax)
legend('Without C.A.','With C. A.','Location','SouthEast')
savefigure([figFileBase '_distances'],'pdf',figDimSmall,flagSaveFigure)
