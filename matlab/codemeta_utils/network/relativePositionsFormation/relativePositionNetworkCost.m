%Compute cost and gradient for relative-positions formation control
function [phi,gradPhi]=relativePositionNetworkCost(E,Tij,TijTruth)
NNodes=max(E(:));
NEdges=size(E,1);
d=size(Tij,1);

phi=0;
gradPhi=zeros(d,NNodes);

for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    phiij=norm(Tij(:,iEdge)-TijTruth(:,iEdge))^2/2;
    gradPhiij=(Tij(:,iEdge)-TijTruth(:,iEdge))*[-1 1];
    phi=phi+phiij;
    gradPhi(:,[iNode jNode])=gradPhi(:,[iNode jNode])+gradPhiij;
end

