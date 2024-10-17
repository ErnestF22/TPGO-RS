function [c,gradc,Hc]=essentialCostNetwork(Ri,Ti,Qij,E)
flagComputeGradient=false;
flagApproximatedHessian=false;
if nargout>1
    flagComputeGradient=true;
    if nargout>2
        flagApproximatedHessian=true;
    end
end

NEdges=size(Qij,3);
c=0;

if flagComputeGradient
    NNodes=size(Ri,3);
    gradc=zeros(6,NNodes);
    nodeIdx=reshape(1:6*NNodes,6,NNodes);
    if flagApproximatedHessian
        Hc=zeros(6*NNodes);
    end
end

for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    
    R1=Ri(:,:,iNode);
    R2=Ri(:,:,jNode);
    T1=Ti(:,iNode);
    T2=Ti(:,jNode);
    Q12=Qij(:,:,iEdge);
    
    if ~flagComputeGradient
        cij=essentialCost(R1,T1,R2,T2,Q12);
    else
        ijNodeIdx=[nodeIdx(:,iNode); nodeIdx(:,jNode)];
        if ~flagApproximatedHessian
            [cij,gradcij]=essentialCost(R1,T1,R2,T2,Q12);
        else
            [cij,gradcij,Hcij]=essentialCost(R1,T1,R2,T2,Q12);
            
            Hc(ijNodeIdx,ijNodeIdx)=Hcij;
        end            
        gradc(ijNodeIdx)=gradc(ijNodeIdx)+gradcij(:);
    end
    c=c+cij;
end
