function CSM15_fRPlot
figDir='../../../papers/network/csm15-tutorial/figures';
b=2;
f=@(x) 1-(1+b.*x)*exp(-b.*x);

fs=setFigFontSize(8);
fn=setFigFont('Times');

funPlot(f,linspace(0,pi))
axis equal
axis([0 3.5 0 1.1])

xlabel('x')
ylabel('f_R(x)')

figFileName=fullfile(figDir,'fR');
savefigure(figFileName,'epsc',[300 150],2)

setFigFontSize(fs);
setFigFont(fn);
