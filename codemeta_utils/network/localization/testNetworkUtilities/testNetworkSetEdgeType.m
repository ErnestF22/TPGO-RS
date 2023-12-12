function t_node=testNetworkSetEdgeType(t_node,EType)
structType=testNetworkDetectStructType(t_node);

switch structType
    case 'single'
        t_node.EType=EType;
    case 'array'
        error('NotImplemented','Symmetrization of edges not implemented for array structure')
end
