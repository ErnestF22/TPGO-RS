function [Eij,QEij]=testNetworkGetRelativeEssential(t_node)
structType=testNetworkDetectStructType(t_node);

switch structType
    case 'single'
        Eij=t_node.Eij;
        QEij=t_node.QEij;
        
    case 'array'
        E=testNetworkGetEdges(t_node);
        NEdges=size(E,1);
        
        Eij=zeros(3,3,NEdges);
        QEij=zeros(6,3,NEdges);
        for iEdge=1:NEdges
            iNode=E(iEdge,1);
            jNode=E(iEdge,2);
            
            Eij(:,:,iEdge)=t_node(iNode).Eij(:,:,jNode);
            QEij(:,:,iEdge)=t_node(iNode).QEij(:,:,jNode);
        end
end
