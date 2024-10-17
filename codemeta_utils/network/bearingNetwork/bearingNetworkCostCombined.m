function [phi,gradPhi]=bearingNetworkCostCombined(...
    EBearings,ERanges,...
    yBearings,yRanges,ygBearings,ygRanges,...
    nyBearings,nyRanges,nygRanges,...
    funsBearings,funsRanges,alpha)

flagComputeGradient=false;
if nargout>1
    flagComputeGradient=true;
end

if ~exist('alpha','var')
    alpha=1;
end

if length(alpha)==1
    alpha=[1 alpha];
end

NNodes=max(EBearings(:));
phiBearing=bearingNetworkCost(EBearings,yBearings,ygBearings,nyBearings,funsBearings);
phiRanges=bearingNetworkCostRanges(ERanges,yRanges,ygRanges,nyRanges,nygRanges,funsRanges,'NNodes',NNodes);
phi=alpha(1)*phiBearing+alpha(2)*phiRanges;

if flagComputeGradient
    gBearing=bearingNetworkCostGradient(EBearings,yBearings,ygBearings,funsBearings);
    gRanges=bearingNetworkCostRangesGradient(ERanges,yRanges,ygRanges,nyRanges,nygRanges,funsRanges,'NNodes',NNodes);
    gradPhi=alpha(1)*gBearing+alpha(2)*gRanges;
end

