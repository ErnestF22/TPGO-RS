function cdc16_comparisonGradientProjection_plots
flagPlotControlEffort=false;
flagPlotConvergence=true;
flagSeparateTrajectories=true;
dataFileBase='cdc16_comparisonGradientProjection_';
allDataFileSets={'standard','oneleader','twoleaders','twoleadersstar'};
%%
close all
figDir='../../../papers/network/cdc16-bearingformation/figures/';
%Column size: 3.25
%Double column size: 7
figDim=[1 3/4]*1.7;
figDimSmall=figDim.*[1 0.75];
flagSaveFigure=2;
ax=[-14 8 -8 9.5];
colors=colorDesaturate([1 0 0; 0 1 0; 0 0 1;0 1 1],0.5);
idxLeaders=[1 2];
fs=setFigFontSize(8);
fn=setFigFont('Times');
%%
for d=2:2
    for iDataFileSet=1:length(allDataFileSets)
        dataFileSet=[num2str(d) 'D_' allDataFileSets{iDataFileSet}];  
        dataFile=[dataFileBase dataFileSet '_data'];
        %%
        figFileBase=fullfile(figDir,[mfilename '_' dataFileSet]);



        %%
        load(dataFile)
        %load('cdc16_comparisonGradientProjection_2D_oneleader_data')
        %load('cdc16_comparisonGradientProjection_2D_twoleaders_data')
        %load('cdc16_comparisonGradientProjection_2D_twoleadersstar_data')
        allMethods=fields(results.x);
        figure(1)
        %fake plots to get legend right
        for iMethod=1:length(allMethods)
            plot([0,1],[0,1],'-','Visible','off','color',colors(iMethod,:))
            hold on
        end
        if ~flagSeparateTrajectories
            for iMethod=1:length(allMethods)
                method=allMethods{iMethod};
                bearingNetworkPlotTrajectories(results.x.(method),'-.','color',colors(iMethod,:))
                hold on
            end
        end
        bearingNetworkPlot(t_node0,'flagPlotCurrent',false)
        hold on
        for iMethod=1:length(allMethods)
            method=allMethods{iMethod};
            bearingNetworkPlot(results.t_node.(method).t_node,'colorCurrent',colors(iMethod,:),...
                'flagPlotDesired',false);
            hold on
            plotPoints(t_node.Titruth(:,idxLeaders(1:NLeaders)),'bs','MarkerSize',8,'MarkerFaceColor',[0 0 1])
        end
        hold off
        axis(ax)
        if d==3
            view(3)
        end
        %legend(allMethods)
        savefigure([figFileBase '_positions'],'pdf',figDim,flagSaveFigure)

        %%
        if flagSeparateTrajectories
            figure(2)
            %fake plots to get legend right
            for iMethod=1:length(allMethods)
                plot([0,1],[0,1],'-','Visible','off','color',colors(iMethod,:))
                hold on
            end
            for iMethod=1:length(allMethods)
                method=allMethods{iMethod};
                bearingNetworkPlotTrajectories(results.x.(method),'-','color',colors(iMethod,:))
                hold on
            end
            hold off
            axis(ax)
            if d==3
                view(3)
            end
            %legend(allMethods)
            savefigure([figFileBase '_trajectories'],'pdf',figDim,flagSaveFigure)
        end        
        %%
        if flagPlotControlEffort
            figure(3)
            for iMethod=1:length(allMethods)
                method=allMethods{iMethod};
                semilogy(results.outputControl.(method).cumulativeNormSq,results.output.(method).dist,...
                    'color',colors(iMethod,:))
                hold on
            end
            hold off
            legend(allMethods)
        end
        %%
        if flagPlotConvergence
            figure(4)
            for iMethod=1:length(allMethods)
                method=allMethods{iMethod};
                semilogy(results.t.(method),results.output.(method).dist,...
                    'color',colors(iMethod,:))
                hold on
            end
            hold off
            %legend(allMethods)
            savefigure([figFileBase '_convergence'],'pdf',figDimSmall,flagSaveFigure)
        end
    end
end
%%

for iColor=1:size(colors,1)
    fprintf('\\definecolor{plot%c}{%.4f,%.4f,%.4f}\n','A'+iColor-1, colors(iColor,:))
end

%%
setFigFontSize(fs);
setFigFont(fn);
