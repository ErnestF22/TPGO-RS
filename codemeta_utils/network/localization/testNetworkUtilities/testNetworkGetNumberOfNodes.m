%function N=testNetworkGetNumberOfNodes(t_node)
%Returns the number of nodes in the structure t_node

function N=testNetworkGetNumberOfNodes(t_node)

switch testNetworkDetectStructType(t_node)
    case 'single'
        N=t_node.NNodes;
    case 'array'
        N=length(t_node);
end
