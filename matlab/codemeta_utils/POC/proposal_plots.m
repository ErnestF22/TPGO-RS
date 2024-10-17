%Plots for NSF RI 18 proposal
%Run POCFactorizationLocalization to generate the data to plot
function proposal_plots
close all
setFigFontSize(7)
setFigFont('Times')

flagDisplayEstimated=false;
figDir='~/BU/proposals/18-NSFRI-LowRankPoseGeometry/figures';
figDim=[3.5 0.8];
flagSaveFile=2;
RPlot=rot_expVec([1.8;0;0]);
load('POCFactorizationLocalization_data.mat')
Ri=multiprod(RPlot,Ri);
Ti=RPlot*Ti;
t_node.Ritruth=Ri;
t_node.Titruth=Ti;

[RProcrustes,RiEstimated]=rotationProcrustes(Ri,RiEstimated,'left');
TiEstimated=RProcrustes*TiEstimated;
TiEstimated=TiEstimated-mean(TiEstimated,2)+mean(Ti,2);
t_node.Ri=RiEstimated;
t_node.Ti=TiEstimated;


figure(1)
    set(gcf, 'Renderer', 'painters')
testNetworkDisplay(t_node,'splitMembers','Ritruth','Titruth',...
    'flagDisplayPoints',false)
if flagDisplayEstimated
    hold on
    testNetworkDisplay(t_node,'splitMembers','Ri','Ti','estimated',...
        'flagDisplayPoints',false)
    hold off
end
xlabel('')
ylabel('')
zlabel('')
savefigure(fullfile(figDir,'network'),'epsc',[figDim(1) 4],flagSaveFile)



figure(2)
semilogy([output.linearResiduals;output.solutionReferenceResiduals]')
ax=axis;
ax(3)=1e-30;
axis(ax)
set(gca,'YTick',[1e-24 1e-12 1])
savefigure(fullfile(figDir,'residuals'),'epsc',figDim,flagSaveFile)

figure(3)
semilogy(output.solutionSvals')
ax=axis;
ax(4)=1e5;
axis(ax)
set(gca,'YTick',[1e-12 1e-6 1])
savefigure(fullfile(figDir,'singular-values'),'epsc',figDim,flagSaveFile)
