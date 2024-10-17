%Returns maximum degree of the network
%function maxDeg=testNetworkGetEdges(t_node)

%%AUTORIGHTS%%

function maxDeg=testNetworkMaxDegree(t_node)
structType=testNetworkDetectStructType(t_node);

switch structType
    case 'single'
        A=t_node.A;
    case 'array'
        A=cat(1,t_node.aij);
end

maxDeg=max(sum(A,1));
