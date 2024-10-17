%Extract arrays of dispersion matrices for the rotation and translations
%function [GammaR,GammaT]=testNetworkGetDispersionMatricesRotationTranslation(t_node)
function [GammaR,GammaT]=testNetworkGetDispersionMatricesRotationTranslation(t_node)
structType=testNetworkDetectStructType(t_node);

switch structType
    case 'single'
        if nargout==1
            GammaR=t_node.dispersionMat;
        else
            GammaR=t_node.dispersionMatR;
            GammaT=t_node.dispersionMatT;
        end
    case 'array'
        E=testNetworkGetEdges(t_node);
        NEdges=size(E,1);
        
        if nargout==1
            Gamma=zeros(6,6,NEdges);
            for iEdge=1:NEdges
                iNode=E(iEdge,1);
                jNode=E(iEdge,2);

                Gamma(:,:,iEdge)=t_node(iNode).dispersionMat(:,:,jNode);
            end
        else
            GammaR=zeros(3,3,NEdges);
            GammaT=zeros(3,3,NEdges);
            for iEdge=1:NEdges
                iNode=E(iEdge,1);
                jNode=E(iEdge,2);

                GammaR(:,:,iEdge)=t_node(iNode).dispersionMatR(:,:,jNode);
                GammaT(:,:,iEdge)=t_node(iNode).dispersionMatT(:,:,jNode);
            end
        end            
end
