%Compute gradient for relative-positions formation control
function gradPhi=relativePositionNetworkCostGradient(E,Tij,TijTruth)
NNodes=max(E(:));
NEdges=size(E,1);
d=size(Tij,1);

gradPhi=zeros(d,NNodes);

for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    gradPhiij=(Tij(:,iEdge)-TijTruth(:,iEdge))*[-1 1];
    gradPhi(:,[iNode jNode])=gradPhi(:,[iNode jNode])+gradPhiij;
end

