function dx=relativePositionNetworkControlEdges(E,Tij,TijTruth)
NEdges=size(E,1);
d=size(Tij,1);

dx=zeros(d,NEdges);

for iEdge=1:NEdges
    dx(:,iEdge)=(Tij(:,iEdge)-TijTruth(:,iEdge));
end
