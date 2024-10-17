function thesis_plotAlbet
close all
flagSaveFigures=true;
figDir='~/Documents/JHU/tron/PhDThesis/figures/';
figFontSize=9;
figDim=[9 5]*savefigureResolution('cm');


oldFontSize=setFigFontSize(figFontSize);

alpha=@(x) x/2*cot(x/2);
beta=@(x) x/2;

albetp=@(x) alpha(x)+sqrt(alpha(x)^2+beta(x)^2);
albetm=@(x) alpha(x)-sqrt(alpha(x)^2+beta(x)^2);

t=linspace(0,pi,200);
plotfun(albetp, t, 'r');
hold on
plotfun(albetm, t, 'b');
hold off
h=legend('\alpha+(\alpha^2+\beta^2)^1^/^2','\alpha-(\alpha^2+\beta^2)^1^/^2', ...
    'Location','SouthWest');
set(h,'Interpreter','TeX')
grid on

setFigFontSize(oldFontSize)
