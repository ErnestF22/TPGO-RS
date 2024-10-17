function journal14_singleIntegrators_plots
%generates plots in the journal paper

fs=setFigFontSize(8);
fn=setFigFont('Times');

for d=[3 2]
    load(['journal14_singleIntegrators_' num2str(d) 'D_standard_data'])

    flagSaveFigure=2;

    figDir='../../../papers/network/TRO14-BearingFormation/figures/';
    figDim=[170,150];
    figDimSmall=[160,100];


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


            figure(5)
            semilogy(t,output.a*180/pi,'b')
            ylabel('Bearing angle error [deg]')
            xlabel('Time')
            switch d
                case 2
                    axis([0 1500 1e-3 1e3])
                    set(gca,'YTick',10.^[-3 -1 1 3])
                case 3
                    axis([0 1500 1e-10 1e0])
            end
            axPos=get(gca,'Position');
            savefigure([figBaseFileName 'angles'],'epsc',figDimSmall,flagSaveFigure)

            figure(9)
            plot(t,abs(output.rd),'b')
            ylabel('Distance error')
            xlabel('Time')
            axis([0 1500 0 12])
            set(gca,'Position',axPos)
            savefigure([figBaseFileName 'rangeDifference'],'epsc',figDimSmall.*[1.055 1],flagSaveFigure)

        end
    end
end

setFigFontSize(fs);
setFigFont(fn);
