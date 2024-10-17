function t_node=testNetworkSetRotTransl(t_node,Ri,Ti)

structType=testNetworkDetectStructType(t_node);

NEdges=size(Ri,3);
Gi=[Ri permute(Ti,[1 3 2]); zeros(1,3,NEdges) ones(1,1,NEdges)];

switch structType
    case 'single'
        t_node.gi=Gi;
    case 'array'
        for iNode=1:length(t_node)
            t_node(iNode).gi=Gi(:,:,iNode);
        end
end
