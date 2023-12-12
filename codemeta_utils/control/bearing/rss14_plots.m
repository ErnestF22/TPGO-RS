function rss14_plots
figDir='../../../papers/network/RSS14-BearingFormation/figures/';
figDim=[300,120];
flagSaveFigure=2;

fs=setFigFontSize(7);

thetac=65/180*pi;
thetas=45/180*pi;

costName='barrier';

cc=cos(thetac);
cs=cos(thetas);
display(cc)
display(cs)

yMax=10;
funs=bearingCostFunctions(costName,cs,cc);

plotfun(funs.f,linspace(cc,1,200))
hold on
plot([cc cc],[0 yMax],'r--')
plot([cs cs],[0 yMax],'r--')
hold off
text(cc,yMax+0.3,'c_c','HorizontalAlignment','center','FontSize',7)
text(cs,yMax+0.3,'c_s','HorizontalAlignment','center','FontSize',7)
axis([0.4 1 0 yMax])

savefigure(fullfile(figDir,'barrier'),'epsc',figDim,flagSaveFigure)
setFigFontSize(fs);
