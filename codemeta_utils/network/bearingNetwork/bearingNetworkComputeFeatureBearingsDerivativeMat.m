function [D,R]=bearingNetworkComputeFeatureBearingsDerivativeMat(E,y,ny,NNodes)
if ~exist('NNodes','var') || isempty(NNodes)
    NNodes=max(E(:,1));
end

d=size(y,1);
NEdges=size(E,1);

%build D
D=kron(diag(ny),eye(d));

%build R
idxNodes=reshape(1:d*NNodes,d,NNodes);
idxEdges=reshape(1:d*NEdges,d,NEdges);
R=zeros(d*NEdges,d*NNodes);
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    
    Py=eye(d)-y(:,iEdge)*y(:,iEdge)';
    R(idxEdges(:,iEdge),idxNodes(:,iNode))=-Py;
end
