%function E=testNetworkGetEdges(t_node)
%Returns the NEdges x 2 matrix E containing the edges extracted from the
%fields t_node().aij

%%AUTORIGHTS%%

function E=testNetworkGetEdges(t_node)
structType=testNetworkDetectStructType(t_node);

switch structType
    case 'single'
        E=t_node.E;
    case 'array'
        A=cat(1,t_node.aij);
        [E1,E2]=find(A~=0);
        E=[E1 E2];
end

