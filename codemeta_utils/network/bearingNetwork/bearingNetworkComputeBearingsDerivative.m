%Compute derivative of bearings for a network
%function [dy,R,B]=bearingNetworkComputeBearingsDerivative(dx,E,y,ny)
%Inputs
%   dx   derivative of agents' locations (i.e., \dot{x})
%   y   current bearing
%   ny  ranges corresponding to the bearings
%   E   network edge list
%Output
%   dy  derivative of y along v
%   R,B matrices to build the differential Dy=inv(R)*B that maps dx to dy
function [dy,R,B]=bearingNetworkComputeBearingsDerivative(dx,E,y,ny)
[d,NNodes]=size(dx);
R=bearingBuildR(ny,d);
B=bearingBuildB(y,NNodes,E);
dy=reshape(R\B*dx(:),d,[]);

function R=bearingBuildR(ny,d)
R=kron(diag(ny),eye(d));

function B=bearingBuildB(y,NNodes,E)
d=size(y,1);
NEdges=size(E,1);
idxNodes=reshape(1:d*NNodes,d,NNodes);
idxEdges=reshape(1:d*NEdges,d,NEdges);
B=zeros(d*NEdges,d*NNodes);
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    
    Py=eye(2)-y(:,iEdge)*y(:,iEdge)';
    B(idxEdges(:,iEdge),idxNodes(:,iNode))=-Py;
    B(idxEdges(:,iEdge),idxNodes(:,jNode))=Py;
end
