function G=rotLocFrobCostMatrix(E,RRel)
NRotations=max(E(:));
NEdges=size(E,1);
I=eye(3);

G=zeros(3*NRotations,3*NRotations);
idxRot=reshape(1:3*NRotations,3,NRotations);
for iEdge=1:NEdges
    idxEdge=[idxRot(:,E(iEdge,1));idxRot(:,E(iEdge,2))];
    G(idxEdge,idxEdge)=G(idxEdge,idxEdge)+rotLocFrobCostPairMatrix(RRel(:,:,iEdge));
end
