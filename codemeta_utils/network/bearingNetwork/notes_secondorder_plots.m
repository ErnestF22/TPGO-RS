function notes_secondorder_plots
figDir='../../../papers/network/secondorderbearingformation/figures';
figDim=[170,150];
figDimSmall=[160,100];

flagSaveFigure=2;
fs=setFigFontSize(8);
fn=setFigFont('Times');

k=[2 2 2];
kNoFeatures=[k(1:2) 0];
allConditions={{'b',k},{'bd',k},...
    {'b',kNoFeatures},{'bd',kNoFeatures}};
    
d=2;
lambda=0;
NConditions=length(allConditions);
for iCondition=1:NConditions
    
    condition=allConditions{iCondition};
    baseFileName=fileNameClean(['notes_secondorder_D' num2str(d) '_L' num2str(lambda)...
        '_' condition{1} '_K' num2str(condition{2},'%d')]);

    disp(baseFileName)
    load([baseFileName 'data'])
    
            %produce figures
            figBaseFileName=fullfile(figDir,baseFileName);
    figure(1)
    bearingNetworkPlot(t_node,'flagPlotFeatureNodes',t_cost.flagUseFeatures)
    hold on
    bearingNetworkPlotTrajectories(x)
    hold off
    savefigure([figBaseFileName 'trajectories'],'epsc',figDim,flagSaveFigure)

    figure(2)
    semilogy(t,output.phi)
    title('Cost')
    savefigure([figBaseFileName 'cost'],'epsc',figDimSmall,flagSaveFigure)

    figure(3)
    plot(t,output.a,'b')
    title('Bearing angle distances')
    savefigure([figBaseFileName 'angles'],'epsc',figDimSmall,flagSaveFigure)

    figure(4)
    plot(t,output.rd,'b')
    title('Range differences')
    savefigure([figBaseFileName 'ranges'],'epsc',figDimSmall,flagSaveFigure)
    
    figure(5)
    plot(t,[squeeze(x(3,:,:)); squeeze(x(4,:,:))],'b-')
    title('Velocities')
    savefigure([figBaseFileName 'velocities'],'epsc',figDimSmall,flagSaveFigure)
end

setFigFontSize(fs);
setFigFont(fn);
    