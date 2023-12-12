function bearingCostGeneralRanges_test
d=2;
N=3;
funsName='squared';

xLandmarks=randn(d,N);
[x,~,~,v]=real_randGeodFun(randn(d,1));
funs=bearingCostFunctions(funsName);

xg=x(0);
[yg,nyg]=bearingCompute(xg,xLandmarks);

check_der(@(t) costAndDer(v,x(t),xLandmarks,yg,nyg,funs),'function',linspace(0,2))

function [c,dc]=costAndDer(v,x,xLandmarks,yg,nyg,funs)
[y,ny]=bearingCompute(x,xLandmarks);

[c,gradc]=bearingCostGeneralRanges(y,yg,ny,nyg,funs);
dc=v'*gradc;
