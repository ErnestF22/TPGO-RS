function [A,b]=translationOnlyLogLikelihoodNetworkMatrix(Ri,Tij,Gamma,E)
NEdges=size(E,1);
NNodes=size(Ri,3);
A=zeros(3*NNodes);
b=zeros(3*NNodes,1);
idxNode=reshape(1:3*NNodes,3,NNodes);
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    
    R1=Ri(:,:,iNode);
    Gamma12=Gamma(:,:,iEdge);
    T12=Tij(:,iEdge);
    
    b1=R1*Gamma12*T12;
    A1=R1*Gamma12*R1';
    
    idxINode=idxNode(:,iNode);
    idxJNode=idxNode(:,jNode);
    
    b(idxINode)=b(idxINode)+b1;
    b(idxJNode)=b(idxJNode)-b1;
    A([idxINode idxJNode],[idxINode idxJNode])=...
        A([idxINode idxJNode],[idxINode idxJNode])+[A1 -A1; -A1 A1];
end
