function gradMat=rotLocFrobCost_gradMatrix(E,R,RRel)
NEdges=size(E,1);
NNodes=size(R,3);
d=size(R,1);
gradMat=zeros(d*NNodes);
idxNodes=reshape(1:d*NNodes,d,NNodes);
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    idxEdge=idxNodes(:,E(iEdge,:));
    
    gradMatPair=rotLocFrobCostPair_gradMatrix(R(:,:,iNode),R(:,:,jNode),RRel(:,:,iEdge));
    
    gradMat(idxEdge,idxEdge)=gradMat(idxEdge,idxEdge)+gradMatPair;
end
