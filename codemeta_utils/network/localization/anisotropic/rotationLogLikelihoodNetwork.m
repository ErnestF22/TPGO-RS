%Cost for localization of rotation using clipped Gaussian at the measurements
%function l=rotationLogLikelihoodNetwork(Ri,Rij,Gammaij,E)
function [l,dl,Hl]=rotationLogLikelihoodNetwork(Ri,Rij,Gammaij,E)
flagComputeGradient=false;
flagApproximatedHessian=false;

if nargout>1
    flagComputeGradient=true;
end
if nargout>2
    flagApproximatedHessian=true;
end

NEdges=size(Rij,3);
l=0;
if flagComputeGradient
    NNodes=size(Ri,3);
    dl=zeros(3,NNodes);
    nodeIdx=reshape(1:3*NNodes,3,NNodes);
    if flagApproximatedHessian
        Hl=zeros(3*NNodes);
    end
end

for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    R1=Ri(:,:,iNode);
    R2=Ri(:,:,jNode);
    R12=Rij(:,:,iEdge);
    Gamma12=Gammaij(:,:,iEdge);
    
    if ~flagComputeGradient
        lij=rotationLogLikelihood(R1,R2,R12,Gamma12);
    else
        ijNodeIdx=[nodeIdx(:,iNode); nodeIdx(:,jNode)];
        if ~flagApproximatedHessian
            [lij,dlij]=rotationLogLikelihood(R1,R2,R12,Gamma12);
        else
            [lij,dlij,Hlij]=rotationLogLikelihood(R1,R2,R12,Gamma12);
            
            Hl(ijNodeIdx,ijNodeIdx)=Hlij;
        end
        dl(ijNodeIdx)=dl(ijNodeIdx)+dlij(:);
    end
    l=l+lij;
end
