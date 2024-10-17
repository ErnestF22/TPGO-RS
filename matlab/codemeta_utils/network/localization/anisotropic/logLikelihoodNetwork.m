function [l,dl,Hl]=logLikelihoodNetwork(Ri,Ti,Rij,Tij,Gammaij,E,EType)
flagComputeGradient=false;
flagApproximatedHessian=false;
if nargout>1
    flagComputeGradient=true;
    if nargout>2
        flagApproximatedHessian=true;
    end
end

NEdges=size(Rij,3);
if ~exist('EType','var')
    EType=ones(NEdges,1);
end

l=0;
if flagComputeGradient
    NNodes=size(Ri,3);
    dl=zeros(6,NNodes);
    nodeIdx=reshape(1:6*NNodes,6,NNodes);
    if flagApproximatedHessian
        Hl=zeros(6*NNodes);
    end
end


for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    
    R1=Ri(:,:,iNode);
    R2=Ri(:,:,jNode);
    T1=Ti(:,iNode);
    T2=Ti(:,jNode);
    R12=Rij(:,:,iEdge);
    T12=Tij(:,iEdge);
    Gamma12=Gammaij(:,:,iEdge);
    EType12=EType(iEdge);
    
    if ~flagComputeGradient
        lij=logLikelihood(R1,R2,T1,T2,R12,T12,Gamma12,EType12); 
    else
        ijNodeIdx=[nodeIdx(:,iNode); nodeIdx(:,jNode)];
        if ~flagApproximatedHessian
            [lij,dlij]=logLikelihood(R1,R2,T1,T2,R12,T12,Gamma12,EType12);
        else
            [lij,dlij,Hlij]=logLikelihood(R1,R2,T1,T2,R12,T12,Gamma12,EType12);
            
            Hl(ijNodeIdx,ijNodeIdx)=Hlij;
        end
            
        dl(ijNodeIdx)=dl(ijNodeIdx)+dlij(:);
    end
    l=l+lij;
end
