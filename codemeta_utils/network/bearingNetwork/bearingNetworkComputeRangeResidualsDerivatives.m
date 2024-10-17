%Compute derivative of range residuals for a network
%function [dq,Bq]=bearingNetworkComputeRangeResidualsDerivatives(dx,E,yg)
%Inputs
%   dx  derivative of agents' locations (i.e., \dot{x})
%   yg  desired bearing
%   E   network edge list
%Output
%   dq  derivative of q along v
%   Bq  differential Dy=B that maps dx to dq
function [dq,Bq]=bearingNetworkComputeRangeResidualsDerivatives(dx,E,yg)
NNodes=size(dx,2);
Bq=bearingBuildBRanges(yg,NNodes,E);
dq=Bq*dx(:);

function Bq=bearingBuildBRanges(yg,NNodes,E)
d=size(yg,1);
NEdges=size(E,1);
idxNodes=reshape(1:d*NNodes,d,NNodes);
Bq=zeros(NEdges,d*NNodes);
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);

    Bq(iEdge,idxNodes(:,iNode))=-yg(:,iEdge);
    Bq(iEdge,idxNodes(:,jNode))=yg(:,iEdge);
end
