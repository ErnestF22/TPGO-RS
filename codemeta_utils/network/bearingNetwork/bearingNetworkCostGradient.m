%Compute gradient of network formation cost
%function gradPhi=bearingNetworkCostGradient(E,y,yg,funs)
function [gradPhi,DgradPhi]=bearingNetworkCostGradient(E,y,yg,funs,ny)
flagComputeDGrad=nargout>1;

NNodes=max(E(:));
NEdges=size(E,1);
d=size(y,1);

gradPhi=zeros(d,NNodes);
if flagComputeDGrad
    idxNodes=reshape(1:d*NNodes,d,NNodes);
    DgradPhi=zeros(d*NNodes,d*NNodes);
end

for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    if ~flagComputeDGrad
        gradPhiij=bearingNetworkCostPair_grad(y(:,iEdge),yg(:,iEdge),funs);
    else
        [gradPhiij,DgradPhiij]=bearingNetworkCostPair_grad(y(:,iEdge),yg(:,iEdge),funs,ny(iEdge));
        DgradPhi(idxNodes(:,E(iEdge,:)),idxNodes(:,E(iEdge,:)))=DgradPhi(idxNodes(:,E(iEdge,:)),idxNodes(:,E(iEdge,:)))+DgradPhiij;
    end
    
    gradPhi(:,[iNode jNode])=gradPhi(:,[iNode jNode])+gradPhiij;
end
