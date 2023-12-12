function cvpr13_plots
close all
%dest='paper';
dest='poster';
flagSaveFile=2;

switch dest
    case 'paper'
        fs=setFigFontSize(7);
        fn=setFigFont('Times');
        figDim=[6 2]*savefigureResolution('cm');
        figDir='~/Documents/UPenn/papers/vision/CVPR13-EssentialManifold/figures/';
    case 'poster'
        fs=setFigFontSize(18);
        fn=setFigFont('Helvetica');
        figDim=2*[6 2]*savefigureResolution('cm');
        figDir='~/Documents/UPenn/presentations/vision/CVPR14-EssentialManifold/figures/';
end

resetRands(1)
Q1=essential_randn();
Q2=essential_randn();

yMin=0;
yMax=14;

f1=@(t) essential_distMinAnglePair_ftFromQ(t,Q1,Q2,'term','first');
f2=@(t) essential_distMinAnglePair_ftFromQ(t,Q1,Q2,'term','second');
f=@(t) essential_distMinAnglePair_ftFromQ(t,Q1,Q2);
f12=@(t) [f1(t); f2(t)];
[tMin,fMin,tBreak1,tBreak2,Q2,tMinAll]=essential_distMinAnglePair(Q1,Q2,1);

t=linspace(-pi,pi,200);
[~,tBreak1]=f1(0);
[~,tBreak2]=f2(0);
plotfun(f,t,'r')
hold on
plotfun(f12,t,'b')
plotfun(f,modAngle(tMinAll),'ro')
plot(modAngle(tBreak1)*[1 1], [yMin yMax], 'k--')
plot(modAngle(tBreak2)*[1 1], [yMin yMax], 'k--')
hold off
axis([-pi pi yMin yMax])

figure(1)
savefigure(fullfile(figDir,'distMinAngle'),'epsc',figDim,flagSaveFile)
setFigFontSize(fs)
setFigFont(fn)
