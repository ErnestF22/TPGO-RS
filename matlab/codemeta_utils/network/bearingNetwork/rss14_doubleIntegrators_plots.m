function rss14_doubleIntegrators_plots
close all
figDir='../../../papers/network/RSS14-BearingFormation/figures';
figDim=[300,200];
figDimSmall=[300,100];

TShort=10;
flagSaveFigure=2;
fs=setFigFontSize(8);
fn=setFigFont('Times');

allCostNamesBearings={'cosine','angleSq'};
lambdas=[0 0.01 0.1 1];
alphaks=[3 6];


NCostBearings=length(allCostNamesBearings);
NLambdas=length(lambdas);
NAlphaks=length(alphaks);
for iCostBearings=1:NCostBearings
    costNameBearings=allCostNamesBearings{iCostBearings};
    disp(['# Cost ' costNameBearings])
    for iLambda=1:NLambdas
        for iAlphak=1:NAlphaks
            alphak=alphaks(iAlphak);
            baseFileName=['bearingNetworkDynamic_' costNameBearings '_'];
            baseFileName=[baseFileName num2str(iLambda) '_' num2str(alphak) '_'];

            figBaseFileName=fullfile(figDir,baseFileName);

            load([baseFileName 'data'],'t','z','x0','x','xg','xFinal','t_node',...
                'phi','funsBearings','c','m','d','funsRanges','q')

            figure(1)
            plot(x0(1,:),x0(2,:),'rx')
            hold on
            plot(squeeze(x(1,:,:))',squeeze(x(2,:,:))','-','color',[1 0.75 0])
            plot([xg(1,:); xFinal(1,:)],[xg(2,:); xFinal(2,:)],'k:')
            bearingNetworkPlot(t_node)
            hold off
            axis([-11 10 -8 8])
            savefigure([figBaseFileName 'trajectories'],'epsc',figDim,flagSaveFigure)

            figure(2)
            semilogy(t,phi)
            savefigure([figBaseFileName 'cost'],'epsc',figDimSmall,flagSaveFigure)

            figure(3)
            semilogy(t,funsBearings.f(c));
            setFigXMax(TShort)
            savefigure([figBaseFileName 'residuals_c'],'epsc',figDim,flagSaveFigure)

            figure(4)
            semilogy(t,funsRanges.f(q));
            setFigXMax(TShort)
            savefigure([figBaseFileName 'residuals_q'],'epsc',figDim,flagSaveFigure)

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

            figure(7)
            plot(t,reshape(x,[],size(x,3)))
            savefigure([figBaseFileName 'positions'],'epsc',figDim,flagSaveFigure)

            figure(8)
            v=z(1:2,:,:);
            plot(t,reshape(v,[],size(v,3)))
            savefigure([figBaseFileName 'velocities'],'epsc',figDim,flagSaveFigure)
        end
    end
end

setFigFontSize(fs);
setFigFont(fn);

function setFigXMax(xMax)
ax=axis();
ax(2)=xMax;
axis(ax)
