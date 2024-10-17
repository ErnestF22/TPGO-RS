function journal14_additional_plots
close all
figDir='../../../papers/network/TRO14-BearingFormation/figures/';
figDim=[300,150];
figDimSmall=[300,100];
flagSaveFigure=2;

fs=setFigFontSize(8);
fn=setFigFont('Times');

switch 1
    case 1
        %moving formation with two leaders
        load('journal14_singleIntegrators_motion_cosine_b_data','t_node','x','idxLeader')
        
        figBaseFileName=fullfile(figDir,'motion_twoLeaders_cosine_');
        figure(1)
        bearingNetworkPlot(t_node,'flagPlotRanges',false,'leader',idxLeader)
        hold on
        bearingNetworkPlotTrajectories(x,'-')  
        plot([t_node.Ti(1,:); t_node.Titruth(1,:)],[t_node.Ti(2,:); t_node.Titruth(2,:)],'k:')        
        hold off
        axis([-43 6 -6.5 19])
        savefigure([figBaseFileName 'trajectories'],'epsc',figDim,flagSaveFigure)
    case 2
        %preconditioner comparison
        figBaseFileName=fullfile(figDir,'preconditioner_cosine_');
        for set={'b','bp'}
            load(['journal14_singleIntegrators_cosine_' set{1} '_data.mat'],'t_node','x')
        
            figure(1)
            bearingNetworkPlot(t_node,'flagPlotRanges',false)
            hold on
            bearingNetworkPlotTrajectories(x)        
            hold off
            axis([-11 10 -8.5 7.5])
            savefigure([figBaseFileName set{1} '_trajectories'],'epsc',figDim,flagSaveFigure)
        end

end

setFigFontSize(fs);
setFigFont(fn);
