%function t_node=testNetworkExecFunOnEdges(t_node,fun,E)
%E must be a NEdges x 2 matrix specifying the edges present in the network.
%If E is not passed, the edges are deduced from the topology defined in the
%fields t_node().aij
%If t_node is of type 'Single':
%Executes a function with prototype t_node=fun(t_node,iEdge), where iEdge
%goes over the edges. The argument E is ignored.
%If t_node is of type 'Array':
%Executes a function with prototype t_node=fun(t_node,iNode,jNode) by
%substituting iNode and jNode with the elements of each row of E.

%%AUTORIGHTS%%

function t_node=testNetworkExecFunOnEdges(t_node,fun,E)

structType=testNetworkDetectStructType(t_node);

switch structType
    case 'single'
        for iEdge=1:t_node.NEdges
            t_node=fun(t_node,iEdge);
        end
    case 'array'
        if ~exist('E','var')
            E=testNetworkGetEdges(t_node);
        end

        NEdges=size(E,1);

        for iEdge=1:NEdges
            iNode=E(iEdge,1);
            jNode=E(iEdge,2);
            t_node=fun(t_node,iNode,jNode);
        end
end


