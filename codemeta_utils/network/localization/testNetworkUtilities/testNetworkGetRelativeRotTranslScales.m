%Extracts arrays relative measured rotation, translation and, optionally, of scales
%function [Rij,tij,lij]=testNetworkGetRelativeRotTranslScales(t_node)
function [Rij,tij,lij]=testNetworkGetRelativeRotTranslScales(t_node)

structType=testNetworkDetectStructType(t_node);

switch structType
    case 'single'
        Gij=t_node.gij;
        if nargout>2
            lij=t_node.lambdaij;
        end
    case 'array'
        E=testNetworkGetEdges(t_node);
        NEdges=size(E,1);
        
        Gij=zeros(4,4,NEdges);
        if nargout>2
            lij=zeros(1,NEdges);
        end
        for iEdge=1:NEdges
            iNode=E(iEdge,1);
            jNode=E(iEdge,2);
            
            Gij(:,:,iEdge)=t_node(iNode).gij(:,:,jNode);
            if nargout>2
                lij(iEdge)=t_node(iNode).lambdaij(jNode);
            end
        end
end

Rij=Gij(1:3,1:3,:);
tij=squeeze(Gij(1:3,4,:));
