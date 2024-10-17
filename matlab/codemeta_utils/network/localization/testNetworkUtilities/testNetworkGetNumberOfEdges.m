%Returns number of edges in the network
function NEdges=testNetworkGetNumberOfEdges(t_node)
structType=testNetworkDetectStructType(t_node);

switch structType
    case 'single'
        NEdges=size(t_node.E,1);
    case 'array'
        A=cat(1,t_node.aij);
        NEdges=sum(A(:)~=0);
end

