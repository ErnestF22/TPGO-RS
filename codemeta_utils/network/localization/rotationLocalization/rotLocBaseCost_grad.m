function gradc=rotLocBaseCost_grad(E,R,RRel,costPair_grad)
NEdges=size(E,1);
gradc=zeros(size(R));
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    gradcPair=costPair_grad(R(:,:,iNode),R(:,:,jNode),RRel(:,:,iEdge));
    gradc(:,:,[iNode jNode])=gradc(:,:,[iNode jNode])+gradcPair;
end
