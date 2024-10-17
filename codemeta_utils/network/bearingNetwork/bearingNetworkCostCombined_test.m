function bearingNetworkCostCombined_test
resetRands()

global funsBearings;
global funsRanges;
costNameBearings='cosine';
costNameRanges='squared';

funsBearings=bearingCostFunctions(costNameBearings);
funsRanges=bearingCostFunctions(costNameRanges);

t_node=bearingNetworkBuildTestNetwork();
[xt,~,x0,v]=real_randGeodFun(t_node.Titruth);
EBearings=t_node.E;
ygBearings=bearingNetworkComputeBearings(x0,EBearings);
ERanges=[1 2;2 1];
[ygRanges,nygRanges]=bearingNetworkComputeBearings(x0,ERanges);

f=@(t) costAndDer(xt(t),EBearings,ERanges,ygBearings,ygRanges,nygRanges,v);
t=linspace(0,1,50);

figure(1)
check_der(f,'function',t)

function [c,dc]=costAndDer(x,EBearings,ERanges,ygBearings,ygRanges,nygRanges,v)
global funsBearings;
global funsRanges;

[yBearings,nyBearings]=bearingNetworkComputeBearings(x,EBearings);
[yRanges,nyRanges]=bearingNetworkComputeBearings(x,ERanges);
[c,gradc]=bearingNetworkCostCombined(...
    EBearings,ERanges,...
    yBearings,yRanges,ygBearings,ygRanges,...
    nyBearings,nyRanges,nygRanges,...
    funsBearings,funsRanges);

dc=v(:)'*gradc(:);
