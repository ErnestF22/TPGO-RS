function journal14_singleIntegrators_motion_plots
%generates plots with moving trajectories in the journal paper

fs=setFigFontSize(8);
fn=setFigFont('Times');

for d=2
    load(['journal14_singleIntegrators_' num2str(d) 'D_motion_data'])
            
    figDir='../../../papers/network/TRO14-BearingFormation/figures/';
    figDim=0.85*[300,200];

    flagSaveFigure=2;

    for iCostBearings=1:NCostBearings
        costNameBearings=allCostNamesBearings{iCostBearings};
        disp(['# Cost ' costNameBearings])

        for iCondition=1:NConditions
            condition=allConditions{iCondition}{1};
            NLeaders=allConditions{iCondition}{2};
            TFinal=allConditions{iCondition}{3};
            baseFileName='journal14_singleIntegrators_';
            if d==3
                baseFileName=[baseFileName '3D_'];
            end
            baseFileName=[baseFileName costNameBearings '_' condition '_'];
            if TFinal~=TLong
                baseFileName=[baseFileName 'T' num2str(TFinal) '_'];
            end
            flagUseRanges=false;
            switch condition
                case 'b'
                    disp('## Bearing formation')
                case 'bp'
                    disp('## Bearing+preconditioner formation')
                case 'bd'
                    disp('## Bearing+distance formation')
                    flagUseRanges=true;
                otherwise
                    error('Test condition for cost not recognized');
            end
            optsPlot={'flagPlotRanges',flagUseRanges};        
            if NLeaders>0
                baseFileName=[baseFileName 'L' num2str(NLeaders) '_'];
            end            

            load([baseFileName 'data'],'t_node','t','x','output','idxLeader')
            optsPlot=[optsPlot 'leader' idxLeader];
            
            

            %produce figures
            figBaseFileName=fullfile(figDir,baseFileName);
            figure(2)
            bearingNetworkPlot(t_node,optsPlot{:})
            hold on
            bearingNetworkPlotTrajectories(x)        
            hold off
            
            sep=2;
            ax=[min(x(1,:))-sep max(x(1,:)+sep) min(x(2,:))-sep max(x(2,:)+sep)];
            axis(ax)
            disp(ax)
            
            
            switch d
                case 2
                    savefigure([figBaseFileName 'trajectories'],'epsc',figDim,flagSaveFigure)
                case 3
                    view(0,90)
                    savefigure([figBaseFileName 'trajectoriesTop'],'epsc',0.5*figDim,flagSaveFigure)
                    view(0,0)
                    savefigure([figBaseFileName 'trajectoriesSide'],'epsc',0.5*figDim,flagSaveFigure)
                    view(-40,6)
                    savefigure([figBaseFileName 'trajectories3Quarters'],'epsc',figDim,flagSaveFigure)
            end

        end
    end
end

setFigFontSize(fs);
setFigFont(fn);