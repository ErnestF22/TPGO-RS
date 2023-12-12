function dx=relativePositionNetworkControl(E,Tij,TijTruth)
NNodes=max(E(:));
NEdges=size(E,1);
d=size(Tij,1);

dx=zeros(d,NNodes);

for iEdge=1:NEdges
    iNode=E(iEdge,1);
    %jNode=E(iEdge,2);
    gradPhiij=(Tij(:,iEdge)-TijTruth(:,iEdge));
    dx(:,iNode)=dx(:,iNode)+gradPhiij;
end

