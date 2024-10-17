function rss14_singleIntegrators_plots
close all
figDir='../../../papers/network/RSS14-BearingFormation/figures';
figDim=[300,150];
figDimSmall=[300,100];

TShort=10;
flagSaveFigure=0;
fs=setFigFontSize(8);
fn=setFigFont('Times');

allCostNamesBearings={'cosine','angleSq'};
NCostBearings=length(allCostNamesBearings);
for iCostBearings=1:1%NCostBearings
    costNameBearings=allCostNamesBearings{iCostBearings};
    disp(['# Cost ' costNameBearings])
    for flagUseRanges=[false true]
        baseFileName=['bearingNetwork_' costNameBearings '_'];
        if ~flagUseRanges
            disp('## Pure bearing formation')
            baseFileName=[baseFileName 'b_'];
        else
            disp('## Bearing+distance formation')
            baseFileName=[baseFileName 'bd_'];
        end

        figBaseFileName=fullfile(figDir,baseFileName);
        
        load([baseFileName 'data'],'t','x0','x','xg','xFinal','t_node',...
            'phi','funsBearings','c','m','d')
        if flagUseRanges
            load([baseFileName 'data'],'funsRanges','q')
        end
            
        
        figure(1)
        plot(x0(1,:),x0(2,:),'rx')
        hold on
        plot(squeeze(x(1,:,:))',squeeze(x(2,:,:))','-','color',[1 0.75 0])
        plot([xg(1,:); xFinal(1,:)],[xg(2,:); xFinal(2,:)],'k:')
        bearingNetworkPlot(t_node,'flagPlotRanges',flagUseRanges)
        hold off
        axis([-11 10 -8.5 7.5])
        savefigure([figBaseFileName 'trajectories'],'epsc',figDim,flagSaveFigure)
        
        figure(2)
        semilogy(t,phi)
        savefigure([figBaseFileName 'cost'],'epsc',figDimSmall,flagSaveFigure)

        figure(3)
        semilogy(t,funsBearings.f(c));
        setFigXMax(TShort)
        savefigure([figBaseFileName 'residuals_c'],'epsc',figDim,flagSaveFigure)
        
        if flagUseRanges
            figure(4)
            semilogy(t,funsRanges.f(q));
            setFigXMax(TShort)
            savefigure([figBaseFileName 'residuals_q'],'epsc',figDim,flagSaveFigure)
        end
        
        figure(5)
        plot(t,m)
        ax=axis();
        ax(3:4)=[-1 1];
        axis(ax)
        savefigure([figBaseFileName 'centroid'],'epsc',figDim,flagSaveFigure)
        
        figure(6)
        plot(t,d)
        setFigXMax(TShort)
        savefigure([figBaseFileName 'distances'],'epsc',figDim,flagSaveFigure)
    end
end

setFigFontSize(fs);
setFigFont(fn);

function setFigXMax(xMax)
ax=axis();
ax(2)=xMax;
axis(ax)
