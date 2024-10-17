function [phi,gradPhi]=bearingNetworkCost(E,y,yg,ny,funs)
NNodes=max(E(:));
NEdges=size(E,1);
d=size(y,1);

phi=0;
gradPhi=zeros(d,NNodes);

for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    [phiij,gradPhiij]=bearingNetworkCostPair(y(:,iEdge),yg(:,iEdge),ny(iEdge),funs);
    phi=phi+phiij;
    gradPhi(:,[iNode jNode])=gradPhi(:,[iNode jNode])+gradPhiij;
end
