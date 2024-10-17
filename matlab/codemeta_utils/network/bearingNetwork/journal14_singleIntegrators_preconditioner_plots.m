function journal14_singleIntegrators_preconditioner_plots
%generates plots for preconditioner in the journal paper

fs=setFigFontSize(8);
fn=setFigFont('Times');

for d=2
    load(['journal14_singleIntegrators_' num2str(d) 'D_preconditioner_data'])
            
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
                optsPlot=[optsPlot 'leader' idxLeader];
            end            

            load([baseFileName 'data'],'t_node','t','x','output')

            %produce figures
            figBaseFileName=fullfile(figDir,baseFileName);
            figure(2)
            bearingNetworkPlot(t_node,optsPlot{:})
            hold on
            bearingNetworkPlotTrajectories(x)        
            hold off
            switch d
                case 2
                    axis([-12 12 -9 7])
                    savefigure([figBaseFileName 'trajectories'],'epsc',figDim,flagSaveFigure)
                case 3
                    axis([-8 8 -8 8 -8 8])
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
