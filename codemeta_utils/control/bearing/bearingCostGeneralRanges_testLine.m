function bearingCostGeneralRanges_testLine

funs=bearingCostFunctions('squared');
xLandmarks=[0;0];
xg=[1;0];

[yg,nyg]=bearingCompute(xg,xLandmarks);
[x,~,~,dx]=real_randGeodFun(xg);

t=linspace(0,2,100);
check_der(@(t) costAndDer(dx,x(t),xLandmarks,yg,nyg,funs),...
    @(t) derB(dx,x(t),xLandmarks,yg,nyg,funs),t);
hold on
plotfun(@(t) funs.f(t*dx'*yg),t,'go')
hold off


function [c,dc]=costAndDer(dx,x,xLandmarks,yg,nyg,funs)
[y,ny]=bearingCompute(x,xLandmarks);
[c,gradc]=bearingCostGeneralRanges(y,yg,ny,nyg,funs);
dc=gradc(:)'*dx(:);

function dc=derB(dx,x,xLandmarks,yg,nyg,funs)
[y,ny]=bearingCompute(x,xLandmarks);
c=bearingComputeCosine(y,yg);
e=ny.*c-nyg;
gradc=-funs.df(e)*yg;
dc=gradc(:)'*dx(:);
