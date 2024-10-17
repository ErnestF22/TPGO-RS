%Compute bearing control field for a network
%function dx=bearingNetworkField(E,y,yg,optField)
function dx=bearingNetworkField(E,y,yg,varargin)
NNodes=max(E(:));
NEdges=size(E,1);
d=size(y,1);

dx=zeros(d,NNodes);

for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    dxij=bearingNetworkFieldPair(y(:,iEdge),yg(:,iEdge),varargin{:});
    dx(:,[iNode jNode])=dx(:,[iNode jNode])+dxij;
end
