%Extract arrays of dispersion matrices for the essential matrices
function Gamma=testNetworkGetDispersionMatriceEssential(t_node)
structType=testNetworkDetectStructType(t_node);

switch structType
    case 'single'
        Gamma=t_node.dispersionMatE;
    case 'array'
        E=testNetworkGetEdges(t_node);
        NEdges=size(E,1);
        
        Gamma=zeros(6,6,NEdges);
        for iEdge=1:NEdges
            iNode=E(iEdge,1);
            jNode=E(iEdge,2);

            Gamma(:,:,iEdge)=t_node(iNode).dispersionMatE(:,:,jNode);
        end
end
